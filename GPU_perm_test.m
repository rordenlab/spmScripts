function GPU_perm_test
%this code demonstrates permutation thresholding when computing millions of correlations
for ii=1:2;
    if (ii==1)
        vv = 70801;
    else
        vv = 239855;
    end
for i=1:20;
    %next - constants for simulations
    kNumRandPerm = 5000; %number of permutations ~5000
    kPcrit = 0.05; %statistical probabiity threshold
    kNumPatients = i *50 ; %in reality this will be >200
    kVoxels = vv; %[53 63 46] for 3mm isotropic and tight bounding box, [182 218 182] = resolution of brain image resliced to 1mm isotropic 70801  238955    1911640
    %next - generate random data
    fprintf('Patients %d  \n',kNumPatients);
    fprintf('Voxels %d  \n',vv );
    les = rand(kNumPatients,kVoxels);
    beh = rand(kNumPatients,1);
    %next - compute and report results
    tic
    [z, threshMin, threshMax] = glm_permSub(les,beh(:,1), kNumRandPerm, kPcrit);
    fprintf('Observed peak %.3f, observed nadir %.3f, Thresholds <%.3f and >%.3f\n',max(z(:)),min(z(:)), threshMin,threshMax);
    toc
end
end
%end perm_test()

function [uncZ, threshMin, threshMax] = glm_permSub(Y, X, nPerms, kPcrit)
%returns uncorrected z-score for all voxels in Y given single column predictor X
% Y: volume data 
% X: single column predictor
%if either X or Y is binomial, results are pooled variance t-test, else correlation coefficient
% nPerms: Number of permutations to compute
% kPcrit: Critical threshold
%Example
% glm_t([1 1 0 0 0 1; 0 0 1 1 1 0]',[1 2 3 4 5 6]') %pooled variance t-test
% 
%inspired by Ged Ridgway's glm_perm_flz
if ~exist('nPerms','var')
    nPerms = 0;
end
if ~exist('kPcrit','var')
    kPcrit = 0.05;
end
% Basics and reusable components
[n f] = size(X); %rows=observations, columns=factors
[nY v] = size(Y); %#ok<NASGU> v is number of tests/voxels
if (f ~= 1), error('glm_permSub is for one column of X at a time (transpose?)'); end; 
if nY ~= n, error('glm_permSub X and Y data sizes are inconsistent'); end
X = [X ones(size(X,1),1)];
c = [1 0]'; %contrast
df = n - rank(X);
pXX = pinv(X)*pinv(X)'; % = pinv(X'*X), which is reusable, because
pX  = pXX * X';         % pinv(P*X) = pinv(X'*P'*P*X)*X'*P' = pXX * (P*X)'
% Original design (identity permutation)
t = glm_quick_t(Y, X, pXX, pX, df, c);
uncZ = spm_t2z(t,df); %report Z scores so DF not relevant
if nPerms < 2
    threshMin = -Inf;
    threshMax = Inf;
    return
end
% Things to track over permutations
peak = nan(nPerms,1);
nadir = nan(nPerms,1);
peak(1) = max(t);
nadir(1) = min(t);
tic
YY= gpuArray(Y);
cc= gpuArray(c); 
GPU_pXX = gpuArray(pXX);
for p = 2:nPerms
    Xp  = X(randperm(n), :);
    GPU_Xp = gpuArray(Xp);
    %pXX = pinv(Xp)*pinv(Xp)'; %??? supposedly not require- reusable?
    pXp = GPU_pXX * GPU_Xp'; % = pinv(Xp)
    tp  = glm_quick_t(YY, GPU_Xp, GPU_pXX, pXp, df, cc);
    c = gather(tp(:));
    peak(p) = max(c);
    nadir(p) = min(c);
end
toc
threshMin = permThreshLowSub (nadir, kPcrit);
threshMax = permThreshHighSub (peak, kPcrit);
threshMin = spm_t2z(threshMin,df); %report Z scores so DF not relevant
threshMax = spm_t2z(threshMax,df); %report Z scores so DF not relevant
%end glm_permSub()


function s = glm_quick_t(y, X, pXX, pX, df, c)
b = pX * y;                     % parameters
y = y - X * b;                  % residuals
s = sum(y .* y);                % sum squared error
s = sqrt(s .* (c'*pXX*c) / df); % standard error
s = c'*b ./ s;                  % t-statistic
%end glm_quick_t()

function thresh = permThreshLowSub (permScores, kPcrit)
permScores = sort(permScores(:));
thresh =permScores(round(numel(permScores) * kPcrit));
%report next most significant score in case of ties
permScores = permScores(permScores < thresh);
if ~isempty(permScores)
    thresh = max(permScores(:));
end
%permThreshLowSub()

function thresh = permThreshHighSub (permScores, kPcrit)
permScores = sort(permScores(:),'descend');
thresh =permScores(round(numel(permScores) * kPcrit));
%report next most significant score in case of ties
permScores = permScores(permScores > thresh);
if ~isempty(permScores)
    thresh = min(permScores(:));
end
%permThreshHighSub()