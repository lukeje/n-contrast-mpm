% array of weighted data structs:
%   weighted.data
%   weighted.fa
%   weighted.TR
%
% Linear fit performed using linearisation assuming short
% TR from Eq (20) of Helms, et al. (2011),
% "Identification of signal bias in the variable flip angle
% method by linear display of the algebraic ernst equation",
% Magn. Reson. Med., https://doi.org/10.1002/mrm.22849

function [A,T1]=weighted2AT1(weightedStructs,relativeB1,mask)

if ~exist('relativeB1','var')
    relativeB1=1;
end

dim = size(weightedStructs(1).data);

if ~exist('mask','var')
    mask = true(dim);
end

Nvoxels=nnz(mask);
Nweighted=numel(weightedStructs);

%% Build design matrix
y=[];
D=[];
for w=1:Nweighted
    weighted = weightedStructs(w);
    
    t=2*tan(relativeB1.*weighted.fa/2);
    
    y=[y,weighted.data(:)./t(:)];                  %#ok<AGROW>
    D=[D,-weighted.data(:).*t(:)/(2*weighted.TR)]; %#ok<AGROW>
end

y=y(mask,:);
D=D(mask,:);

%% Solve for R2s and coefficients
a=zeros(Nvoxels,1);
t1=zeros(Nvoxels,1);
parfor n=1:Nvoxels
    b=y(n,:)/[ones(1,Nweighted);D(n,:)];
    a(n)=b(1);
    t1(n)=b(2);
end

%% Output
A=nan(dim);
A(mask)=a;

T1=nan(dim);
T1(mask)=t1;

end
