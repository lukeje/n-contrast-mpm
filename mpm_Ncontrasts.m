function mpm_Ncontrasts(contrasts, b1map, outdir, threshold)
% Compute T1, PD and R2* from any number of flip angle acquisitions
%
% Reads metadata from sidecar json files
%
% Requires SPM and the hMRI toolbox
%
% Data must have been registered and resliced to the same space outside this script!
%
% contrasts: a cell array of string arrays giving paths to data, e.g.
%     contrasts = {["PDw_e1.nii", "PDw_e2.nii", "PDw_e3.nii"], 
%                  ["T1w_e1.nii", "T1w_e2.nii", "T1w_e3.nii"],
%                  ["ern_e1.nii", "ern_e2.nii", "ern_e3.nii"]};
% outdir:    output directory
% b1map:     B1map resliced to MPM space
% threshold: threshold to mask low intensity data for T1 map calculation

Vref = spm_vol(char(contrasts{1}(1)));

%% Get files
% Extract metadata
for c = 1:length(contrasts)
    currentFiles = contrasts{c}(:);
    nEchoes = length(currentFiles);
    weightedData(c).TE = zeros([1,nEchoes]); %#ok<AGROW>
    weightedData(c).TR = zeros([1,nEchoes]); %#ok<AGROW>
    weightedData(c).fa = zeros([1,nEchoes]); %#ok<AGROW>
    for echo = 1:nEchoes
        imageFile = fullfile(currentFiles(echo).folder,currentFiles(echo).name);
        fid = fopen(strrep(imageFile,"nii","json"), 'r');
        bidsJson = jsondecode(fscanf(fid,"%s"));
        bidsJson = bidsJson.acqpar;
        fclose(fid);
        weightedData(c).TE(echo) = bidsJson.EchoTime*1e-3;
        weightedData(c).TR(echo) = bidsJson.RepetitionTime*1e-3;
        weightedData(c).fa(echo) = deg2rad(bidsJson.FlipAngle);
        assert(weightedData(c).TR(echo)==weightedData(c).TR(1), "TR must match within a contrast!")
        assert(weightedData(c).fa(echo)==weightedData(c).fa(1), "Flip angle must match within a contrast!")
    end
end

%% Fit R2*
R2star = nan(Vref.dim);
extrapolated = cell(length(contrasts),1);
for c = 1:length(contrasts)
    extrapolated{c} = nan(Vref.dim);
end
for z=1:Vref.dim(3) % process data by slice
    for c = 1:length(contrasts)
        currentFiles = char(contrasts{c}(:)); % SPM requires char array
        nEchoes = length(currentFiles);
        weightedData(c).data = nan([Vref.dim(1:2),nEchoes]);
        for echo = 1:nEchoes
            imageFile = fullfile(currentFiles(echo).folder,currentFiles(echo).name);
            V=spm_vol(imageFile);
            weightedData(c).data(:,:,echo) = spm_slice_vol(V,spm_matrix([0 0 z]),Vref.dim(1:2),0);
        end
        weightedData(c).data(weightedData(c).data<eps) = eps;
    end
    [R2star(:,:,z),extrapolatedz] = hmri_calc_R2s(weightedData,"WLS1");
    for c = 1:length(contrasts)
        extrapolated{c}(:,:,z) = extrapolatedz{c};
    end
end

Vout = Vref;
Vout.dt(1) = spm_type('float32');
Vout.fname = char(fullfile(outdir,"R2star.nii"));
spm_write_vol(Vout,R2star);

TEzerofile = cell(length(contrasts),1);
for c = 1:length(contrasts)
    TEzerofile{c} = fullfile(char(outdir),spm_file(spm_file(char(contrasts{c}(1)),'suffix',['con-',num2str(c),'_TEzero']),'filename'));
    Vout.fname = TEzerofile{c};
    spm_write_vol(Vout,extrapolated{c});
end

clear R2star extrapolated

%%
for c = 1:length(contrasts)
    dat0(c).TR = weightedData(c).TR(1); %#ok<AGROW>
    dat0(c).fa = weightedData(c).fa(1); %#ok<AGROW>
end

clear weightedData

%%
T1 = nan(Vref.dim);
A  = nan(Vref.dim);
for z = 1:Vref.dim(3) % process data by slice
    B1 = hmri_read_vols(spm_vol(char(fullfile(b1map.folder,b1map.name))),Vref,z,3)*0.01;
    for c = 1:length(contrasts)
        dat0(c).data = hmri_read_vols(spm_vol(TEzerofile{c}),Vref,z,3);
    end

    mask = t1w0.data>threshold;

    [A(:,:,z),T1(:,:,z)]=weighted2AT1(dat0,B1,mask);
end

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