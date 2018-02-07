function fix_qform (fnm)
%set qform to zero https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
% fnm : nifti image to fix
%Examples
% fix_qform; %use graphical interface
% fix_qform('test.nii');

if ~exist('fnm','var')  %fnm not specified
   [A,Apth] = uigetfile({'*.nii;';'*.*'},'Select image to patch');
   fnm = fullfile(Apth, A);
end;
[pth,nam,ext] = fileparts( fnm);
if ~strcmpi(ext,'.nii'), error('Requires uncompressed .nii files'); end;
qOffsetBytes = 252;
%check that file needs fixing
fid = fopen(fnm);
fseek(fid,qOffsetBytes,'bof');
qform_code = fread(fid,1,'int16');
fclose(fid);
if qform_code == 0, error('qform already set to zero %s', fnm); end;
%create copy of image and modify
outnm = fullfile(pth,['z',nam,ext]);
copyfile(fnm,outnm);
fid = fopen(fnm);
[data,count]=fread(fid, 'uint8');
fclose(fid);
fid = fopen(outnm,'w');
%modify both bytes of 16-bit qform_code
data(qOffsetBytes) = 0;
data(qOffsetBytes+1) = 0;
fwrite(fid,data);
fclose(fid);
%end fix_qform()
