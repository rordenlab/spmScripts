function [tr] =  nii_getTR(fmriname);
%Returns Repeat Time in seconds often stored in raw NIfTI header (though removed during processing)
%  fmriname : name of NIFTI format fMRI image
%Example
% nii_getTR('fmri.nii');

if nargin <1 %no input: select fMRI file[s]
 fmriname = spm_select(1,'image','Select fMRI volume');
end
[pth,nam,ext,vol] = spm_fileparts( deblank (fmriname(1,:)));
fmrinameVol1 = fullfile(pth,[ nam, ext,',1']); 
hdr = spm_vol(fmrinameVol1);
exist('hdr.private.timing.tspacex','var')
if isfield(hdr(1,1).private.timing,'tspace')
  tr = hdr(1,1).private.timing.tspace;
else
  fprintf('%s error: unable to determine TR for image %s\n',mfilename,fmriname); 
end

