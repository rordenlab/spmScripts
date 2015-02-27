function nii_nii2voi (V);
% Converts a NIfTI image to a MRIcron VOI file.
%  A VOI file is a GZipped NIfTI image with the extension '.VOI'
%  MRIcron uses VOI files for binary volumes of interest (e.g. lesion maps)
%  Currently, only works for .nii images, not .hdr/.img images
% Example
%   nii_nii2voi('C:\irate\chrisr.nii');

if nargin <1 %no files
 V = spm_select(inf,'image','Select images to delete');
end;
for i=1:size(V,1)
  ref = deblank(V(i,:));
  [pth,nam,ext] = spm_fileparts(ref);
  if (exist(ref) ~= 2) 
 	fprintf('nii_nii2voi warning unable to find file %s.\n',ref);
	return;  
  end;
  upext = upper(ext);
  if (strcmp(upext,'.IMG') ||  strcmp(upext,'.HDR'))
 	fprintf('nii_nii2voi warning: currently only supports .nii files, not .hdr/.img pairs');
	return; 
  end;
  gzip(ref);
  movefile(fullfile(pth,[nam ext '.gz']), fullfile(pth,[nam '.voi']));
end; %for each file

