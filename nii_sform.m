function nii_sform (pth)
%Report NIfTI header for .nii in path
% pth : folder to convert
%Examples
% nii_sform(pwd)
%To convert images
% dicm2nii(fullfile(pwd,'In'), fullfile(pwd,'Out'),0);

if ~exist('pth', 'var'), pth = pwd; end;
f = dir(fullfile(pth,'*.nii'));
for i = 1: numel(f) 
    fnm = fullfile(pth, f(i).name);
    h = spm_vol(fnm);%
    m = h(1).mat;
    [~, fnm] = fileparts(fnm);
    fprintf('%s = [\n %g %g %g %g;\n %g %g %g %g;\n %g %g %g %g;\n 0 0 0 1]\n', fnm, ...
        m(1,1), m(1,2), m(1,3), m(1,4), ...
        m(2,1), m(2,2), m(2,3), m(2,4), ...
        m(3,1), m(3,2), m(3,3), m(3,4)); 
end