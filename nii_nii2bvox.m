function nii_nii2bvox(filename)
%convert NIfTI-format volume of voxels to a Blender voxel-format
% filename   : name of NIfTI image to convert
%See
% http://pythology.blogspot.com/2014/08/you-can-do-cool-stuff-with-manual.html
% smooth     : blur image
%Examples
% nii_nii2bvox %use graphical interface
% nii_nii2bvox('mni152_2009_256.nii.gz');
%Chris Rorden 2015, BSD license
% https://github.com/bonilhamusclab/MRIcroS/blob/master/%2BfileUtils/readVox.m

if ~exist('filename', 'var')  
   [A,Apth] = uigetfile({'*.nii;*.hdr;*.gz;'},'Select NIFTI image', 'MultiSelect', 'off');
   filename = strcat(Apth,char(A));
end;
if ~exist(filename,'file'), error('Unable to find image %s', filename); end;
[p,n,x] = fileparts(filename);
if strcmpi(x,'.gz')
    filename = char(gunzip(filename));
end
hdr = spm_vol(filename);
img = spm_read_vols(hdr);
if strcmpi(x,'.gz')
    delete(filename ); %FSL hates filename.nii co-existing with filename.nii.gz
end
%scale intensity range
img(isnan(img)) = 0; %remove not-a-numbers
img = img - min(img(:)); %translate so minumum value is 0
img = img / max(img(:)); %scale so image range now 0..1
%write data
outname = fullfile(p, [n '.bvox']);  
fid = fopen(outname,'wb');
fwrite(fid, size(img,1), 'uint32'); %X dimension
fwrite(fid, size(img,2), 'uint32'); %Y dimension
fwrite(fid, size(img,3), 'uint32'); %Z dimension
fwrite(fid, 1, 'uint32'); %one volume
img = single(img(:));
fwrite(fid, img, 'single');
fclose(fid);
%nii_nii2bvox()
