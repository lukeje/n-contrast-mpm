function mpm_Ncontrasts(contrasts, b1map, outdir, threshold, r2sfamethod)
% Compute T1, PD and R2* from any number of flip angle acquisitions
%
% Reads metadata from sidecar json files
%
% Requires SPM and the hMRI toolbox
%
% Data must have been registered and resliced to the same space outside this script!
%
% contrasts:   a cell array of string arrays giving paths to data, e.g.
%     contrasts = {["PDw_e1.nii", "PDw_e2.nii", "PDw_e3.nii"], 
%                  ["T1w_e1.nii", "T1w_e2.nii", "T1w_e3.nii"],
%                  ["ern_e1.nii", "ern_e2.nii", "ern_e3.nii"]};
% outdir:      output directory
% b1map:       path(s) of B1map resliced to MPM space; can be one per contrast or just one
% threshold:   threshold to mask low intensity data for T1 map calculation
% r2sfamethod: method to account for flip-angle dependence of R2* estimation
%              ('none' or 'linear').

Vref = spm_vol(char(contrasts{1}(1)));
Vref = Vref(1); % in case the first file is 4D

% ensure there is a B1 map for each contrast
if (ischar(b1map) && size(b1map,1)==1) || isscalar(b1map)
    b1map = repmat(b1map,length(contrasts),1);
end

%% Get files
% Extract metadata
for c = length(contrasts):-1:1 % allocate backwards to get rid of matlab warnings about needing preallocation
    currentFiles = contrasts{c}(:);
    nFiles = length(currentFiles);
    V{c} = spm_vol(cellstr(currentFiles));
    V{c} = cat(1,V{c}{:});
    nEchoes = length(V{c});
    weightedData(c).TE = [];
    weightedData(c).TR = [];
    fa = [];
    for fIdx = nFiles:-1:1
        imageFile = currentFiles(fIdx);
        fid = fopen(strrep(imageFile,"nii","json"), 'r');
        bidsJson = jsondecode(fscanf(fid,"%s"));
        fclose(fid);

        if isfield(bidsJson,'acqpar') % hMRI toolbox style
            bidsJson = bidsJson.acqpar;
            weightedData(c).TE(:,fIdx) = bidsJson.EchoTime*1e-3;
            weightedData(c).TR(fIdx) = bidsJson.RepetitionTime*1e-3;
        else
            weightedData(c).TE(:,fIdx) = bidsJson.EchoTime;
            if isfield(bidsJson,'RepetitionTimeExcitation') % new BIDS
                weightedData(c).TR(fIdx) = bidsJson.RepetitionTimeExcitation;
            else % old BIDS
                weightedData(c).TR(fIdx) = bidsJson.RepetitionTime;
            end
        end
        fa(fIdx) = deg2rad(bidsJson.FlipAngle);
        
        assert(weightedData(c).TR(fIdx) == weightedData(c).TR(end), "TR must match within a contrast!")
        assert(fa(fIdx) == fa(end), "Flip angle must match within a contrast!")
    end
    weightedData(c).TE = weightedData(c).TE(:);
    assert(length(weightedData(c).TE) == nEchoes, "There must be as many TEs as files per contrast!")
    weightedData(c).TR = weightedData(c).TR(end);
    weightedData(c).fanom = fa(end);
end

%% Fit R2*
R2star = nan(Vref.dim);
DeltaR2star = nan(Vref.dim);
extrapolated = cell(length(V),1);
for c = 1:length(V)
    extrapolated{c} = nan(Vref.dim);
end
spm_progress_bar('Init',Vref.dim(3),'R2* fit');
for z=1:Vref.dim(3) % process data by slice
    for c = 1:length(V)
        nEchoes = length(V{c});
        weightedData(c).data = nan([Vref.dim(1:2),nEchoes]);
        for fIdx = 1:nEchoes
            weightedData(c).data(:,:,fIdx) = spm_slice_vol(V{c}(fIdx),spm_matrix([0 0 z]),Vref.dim(1:2),0);
        end
        weightedData(c).data(weightedData(c).data<eps) = eps;

        B1 = hmri_read_vols(spm_vol(char(b1map(c,:))),Vref,z,3)*0.01;
        weightedData(c).fa = weightedData(c).fanom*B1;
    end

    [R2star(:,:,z),extrapolatedz,DeltaR2star(:,:,z)] = weighted2R2s(weightedData,"WLS1",r2sfamethod);

    for c = 1:length(contrasts)
        extrapolated{c}(:,:,z) = extrapolatedz{c};
    end
    spm_progress_bar('Set',z);
end
spm_progress_bar('Clear');

Vout = Vref;
Vout.dt(1) = spm_type('float32');
Vout.fname = char(fullfile(outdir,"R2starMap.nii"));
spm_write_vol(Vout,R2star);

Vout.fname = char(fullfile(outdir,"DeltaR2starMap.nii"));
spm_write_vol(Vout,DeltaR2star);

TEzerofile = cell(length(contrasts),1);
for c = 1:length(contrasts)
    file = char(contrasts{c}(1));
    TEzerofile{c} = fullfile(char(outdir),spm_file(spm_file(file,'suffix',['_con-',num2str(c),'_TEzero']),'filename'));
    Vout.fname = TEzerofile{c};
    spm_write_vol(Vout,extrapolated{c});

    % also output a sidecar file
    json = jsondecode(fileread(strrep(file,".nii",".json")));
    json.EchoTime = 0;
    json.EchoNumber = 0;
    fid = fopen(strrep(TEzerofile{c},".nii",".json"), 'w');
    fprintf(fid, "%s", jsonencode(json, "PrettyPrint",true));
    fclose(fid);
end

clear R2star extrapolated

%%
for c = length(contrasts):-1:1
    dat0(c).TR    = weightedData(c).TR;
    dat0(c).fanom = weightedData(c).fanom;
end

%%
T1 = nan(Vref.dim);
A  = nan(Vref.dim);
spm_progress_bar('Init',Vref.dim(3),'T1 fit');
for z = 1:Vref.dim(3) % process data by slice
    for c = 1:length(contrasts)
        dat0(c).data = hmri_read_vols(spm_vol(TEzerofile{c}),Vref,z,3);

        B1 = hmri_read_vols(spm_vol(char(b1map(c,:))),Vref,z,3)*0.01;
        dat0(c).fa = B1*dat0(c).fanom;
    end

    mask = dat0(1).data>threshold;

    [A(:,:,z),T1(:,:,z)]=weighted2AT1(dat0,1,mask);

    spm_progress_bar('Set',z);
end
spm_progress_bar('Clear');

VT1 = Vout;
VT1.fname = char(fullfile(outdir,"T1map.nii"));
spm_write_vol(VT1,T1);

VR1 = Vout;
VR1.fname = char(fullfile(outdir,"R1map.nii"));
spm_write_vol(VR1,1./T1);

VA = Vout;
VA.fname = char(fullfile(outdir,"Amap.nii"));
spm_write_vol(VA,A);

end