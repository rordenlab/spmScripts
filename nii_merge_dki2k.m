function nii_merge_dki2k(pth)
%concatenate ["DTI2K_"+"DTI_"] and ["DTI2Krev_"+"DTIrev_"]
% pth : folder with images
%Examples
% nii_merge_dti2k %gui
% nii_merge_dti2k(pwd)

if ~exist('pth','var')
    pth = uigetdir(pwd);
end;

DTI = findSub(pth, 'DTI_*.nii');
DTI2K = findSub(pth, 'DTI2K_*.nii');
DTIrev = findSub(pth, 'DTIrev_*.nii');
DTI2Krev = findSub(pth, 'DTI2Krev_*.nii');
nii_merge_dki({DTI, DTI2K},false, false);
nii_merge_dki({DTIrev, DTI2Krev},false, true);

function fnm = findSub(pth, fnm)
nam = fullfile(pth, fnm);
d = dir(nam);
if isempty(d),
    nam = fullfile(pth, fnm, '.gz');
    d = dir(nam);
end;

if isempty(d), error('Unable to find .nii/.nii.gz named %s', nam); end;
if numel(d) > 1, error('Unable to distinguish between multiple files %s', nam); end;
fnm = fullfile(pth, d.name);
%end findSub()