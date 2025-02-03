function [R2s,extrapolated,DeltaR2s]=weighted2R2s(weighted_data,fitmethod,famethod)
% R2* estimation using an implementation of the ESTATICS
% model (Weiskopf2014). Can utilise weighted least squares (WLS) instead of
% the original ordinary least squares (OLS; Weiskopf2014) to account
% for the heteroscedasticity of log transformed data (Edwards2022).
% Can also estimate a linear flip angle dependence of R2* (Milotta2023).
%
% Input:
%   array of structures (one per contrast) in the form:
%     weighted(contrast).data (NvoxelsX x NvoxelsY x ... x Nechoes)
%     weighted(contrast).TE  (1 x Nechoes)
%   -Voxels must correspond between the weightings (i.e. the images should
%    have been resliced to the same space), but the sampled TEs may be
%    different.
%   -Nechoes must be at least 2 for each weighting.
%   -Because log(0) is ill-defined, zero values in any voxel will result
%    in NaN output for that voxel. To avoid potentially biasing the data,
%    we do not modify the input in any way to avoid this, and leave it to
%    the user to decide how to handle this case, e.g. by removing the
%    corresponding voxels from the input data or replacing zeroes with a
%    small positive number.
%   -For famethod 'linear', the structures should also have an 'fa' field
%    giving the local flip angle in radians or degrees (output will be in
%    reciprocal units) and fanom giving nominal flip angle for the initial
%    OLS fit to get initial weights.
%
%   fitmethod:
%     string stating which fitting method to use.
%     -Options are: 'WLS[N]' (log-linear weighted least squares estimate
%                             with '[N]' iterations, where '[N]' is 1, 2,
%                             or 3; uses OLS signal estimates for initial
%                             weights)
%     -The weights in the WLS case depend on the unknown true signal
%      intensities. These weights can be iteratively updated using the
%      estimated signal intensities. However the benefit of
%      iteratively updating the weights has been found to be relatively
%      small for typical MPM data, and so 'wls1' seems to be sufficient
%      to improve R2* map quality over OLS (Edwards2022).
%
%   famethod:
%     string stating whether to take into account the flip-angle
%     dependence of R2*.
%     -Options are: 'none'   (Weiskopf2014)
%                   'linear' (Milotta2023)
%
% Outputs:
%   R2s (NvoxelsX x NvoxelsY x ...): the voxelwise-estimated
%       common R2* of the weightings.
%   extrapolated: cell array containing data extrapolated to TE=0 in the
%       same order as the input (e.g. matching contrast order).
%   DeltaR2s (NvoxelsX x NvoxelsY x ...): the flip-angle dependent
%       common R2* of the weightings (zero if famethod='none')
%
% Examples:
%   WLS estimate of R2* from PD-weighted data with 3 iterations:
%       R2s = hmri_calc_R2s(struct('data',PDw,'TE',PDwTE),'WLS3');
%
%   WLS-ESTATICS estimate of common R2* of PD, T1, and MT-weighted data
%   with 1 iteration:
%       R2s = hmri_calc_R2s([struct('data',PDw,'TE',PDwTE),...
%           struct('data',T1w,'TE',T1wTE), struct('data',MTw,'TE',MTwTE)],'WLS1');
%
% References:
%   Weiskopf et al. Front. Neurosci. (2014), "Estimating the apparent
%     transverse relaxation time (R2*) from images with different contrasts
%     (ESTATICS) reduces motion artifacts",
%     https://doi.org/10.3389/fnins.2014.00278
%   Edwards et al. Proc. Int. Soc. Magn. Reson. Med. (2022), "Robust and
%     efficient R2* estimation in human brain using log-linear weighted
%     least squares"
%   Milotta et al. Magn. Reson. Med. (2023), "Mitigating the impact of
%     flip angle and orientation dependence in single compartment R2*
%     estimates via 2-pool modeling."
%     https://doi.org/10.1002/mrm.29428

assert(isstruct(weighted_data),'hmri:structError',['inputs must be structs; see help ' mfilename])

dims=size(weighted_data(1).data);
Nvoxels=prod(dims(1:end-1));
Nweighted=numel(weighted_data);

% default to classic ESTATICS
if ~exist('famethod','var') || isempty(famethod), famethod='none'; end

% Account for extra fitting parameter
switch lower(famethod)
    case 'none'
        wBegin=2;
    case 'linear'
        wBegin=3;
    otherwise
        error("flip angle-dependence model '%s' unknown", famethod)
end

%% Build regression arrays
% Build design matrix
D=[];
fa=[];
for w=1:Nweighted
    d=zeros(length(weighted_data(w).TE),Nweighted+wBegin-1);
    d(:,1)=-weighted_data(w).TE;
    d(:,wBegin+w-1)=1;
    switch lower(famethod)
        case 'linear'
            assert(isfield(weighted_data(w),'fa'),'flip angle must be present in weighted_data.fa for each struct!') 
            assert(isfield(weighted_data(w),'fanom'),'nominal flip angle must be present in weighted_data.fanom for each struct!')
            assert(isscalar(weighted_data(w).fanom),'please specify weighted_data.fanom as a scalar, i.e. the nominal flip angle')
            d(:,2)=d(:,1)*weighted_data(w).fanom; % nominal FA used for fast OLS initialisation
            fa = [fa,weighted_data(w).fa(:)]; % save flip angles as array which can be sliced for parfor
    end
    D=[D;d]; %#ok<AGROW>
end

% Build response variable vector
y=[];
for w=1:Nweighted
    
    nTEs=length(weighted_data(w).TE);
    assert(nTEs>1,'each weighting must have more than one TE')
    
    localDims=size(weighted_data(w).data);
    assert(localDims(end)==nTEs,'echoes must be in the final dimension')
    assert(prod(localDims(1:end-1))==Nvoxels,'all input data must have the same number of voxels');
    
    rData=reshape(weighted_data(w).data,Nvoxels,nTEs);
    
    % log(0) is not defined, so warn the user about zeroes in their data
    % for fitting methods involving a log transform.
    % The warning can be disabled with "warning('off','hmri:zerosInInput')"
    if any(rData(:)==0)&&~contains(lower(fitmethod),'nlls')
        warning('hmri:zerosInInput',[...
            'Zero values detected in some voxels in the input data. This ',...
            'will cause estimation to fail in these voxels due to the log ',...
            'transform. If these voxels are background voxels, consider ',...
            'removing them from the input data matrices. ',...
            'Zero values which occur only at high TE in voxels of interest ',...
            'could be replaced with a small positive number, e.g. eps(1) ',...
            '(if the data magnitudes are ~1) or 1 if the data are ',...
            'integer-valued. Note: Care must be taken when replacing ',...
            'values, as this could bias the R2* estimation.']);
    end
    
    y=[y;rData.']; %#ok<AGROW>
end

%% Estimate R2*
switch lower(fitmethod)
    case {'wls1','wls2','wls3'}
        % Number of WLS iterations is specified using 'WLS[N]', where '[N]'
        % is a positive integer
        r=regexp(lower(fitmethod),'^wls(\d+)$','tokens');
        niter=str2double(r{1}{1});
        
        logy=log(y);
        
        % Use OLS estimates for initial weights
        y0=exp(D*OLS(logy,D));
        
        % Loop over voxels
        parfor n=1:size(y,2)
            Dloc = D;
            if strcmp(famethod,'linear')
                for w=1:Nweighted
                    Dloc(2,Dloc(wBegin+w-1)==1) = fa(n,w);
                end
            end
            beta(:,n)=WLS(logy(:,n),Dloc,y0(:,n),niter);
        end
        beta(wBegin:end,:)=exp(beta(wBegin:end,:));
    otherwise
        error("fitting method '%s' not recognised", fitmethod)
end

%% Output
% extra unity in reshape argument avoids problems if size(dims)==2.
R2s=reshape(beta(1,:),[dims(1:end-1),1]);

switch lower(famethod)
    case 'none'
        DeltaR2s=zeros([dims(1:end-1),1]);
    case 'linear'
        DeltaR2s=reshape(beta(2,:),[dims(1:end-1),1]);
end

% Extrapolate weightings to TE=0
if nargout>1
    extrapolated=cell(size(weighted_data)); % cell element per contrast
    for w=1:Nweighted
        extrapolatedData=beta(wBegin+w-1,:);
        extrapolated{w}=reshape(extrapolatedData(:),[dims(1:end-1),1]);
    end
end

end

%% Fitting methods
function beta=OLS(y,D)
% allows for vectorized voxel processing

beta=(D'*D)\(D'*y);

end

function beta=WLS(y,D,y0,niter)
% y0 estimate needed for initial weights, could be raw y values or estimate
% from OLS.
% This function only allows single voxel processing!

% Fix number of iterations to avoid checking convergence of each voxel
for m=1:niter
    % weights are updated using latest parameter estimates
    W=diag(y0.*conj(y0));
    W=W./trace(W); % normalisation
    
    % WLS estimate
    beta=(D'*W*D)\(D'*W*y);
    
    y0=exp(D*beta);
end

end
