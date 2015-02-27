function nii_delete (V);
% This script deletes NIfTI images
%  If passed .nii then deletes one file,
%  If passed .hdr or .img it deletes BOTH files
% Example
%   nii_midline('C:\irate\chrisr.nii');

if nargin <1 %no files
 V = spm_select(inf,'image','Select images to delete');
end;

for i=1:size(V,1)
  ref = deblank(V(i,:));
  [pth,nam,ext] = spm_fileparts(ref);
  fname = fullfile(pth,[ nam ext]); %the T1 image has no prefix
  if (exist(fname) ~= 2) 
 	fprintf('nii_delete warning unable to find file %s.\n',fname);
	return;  
  end;
  delete(fname);

  upext = upper(ext);
  if strcmp(upext,'.IMG')
 	fname = fullfile(pth,[ nam '.hdr']);
 	if (exist(fname) == 2)
 	  delete(fname);
 	end; 
  end;

  if strcmp(upext,'.HDR')
 	fname = fullfile(pth,[ nam '.img']);
 	if (exist(fname) == 2)
 	  delete(fname);
 	end; 
  end;
end; %for each file in V

