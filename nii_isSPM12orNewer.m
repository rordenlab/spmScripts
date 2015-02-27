function nii_isSPM12orNewer
%check that SPM is installed and is at least release 6225
if exist('spm','file') ~= 2, error('Please install SPM12 or later'); end;
[v,r] = spm('Ver','',1); r = str2double(r); %#ok<ASGLU>
if r < 6225, error('Please update your copy of SPM to recent SPM12 or later'); end;
%end isSPM12orNewer()