function nii_zeronan(donor, recipient)
%Voxels with NaN's are replaced with zero, output has 'z' prefix
% donor : image that donates a header
% recipient : image that receives the header
%Example
% nii_swaphdr('oT2_LM1052.nii','LS_LM1052.nii');

if ~exist('donor','var')
	donor = spm_select(1,'image','Select donor image'); 
end
if ~exist('recipient','var')
	recipient = spm_select(1,'image','Select recipient image'); 
end

dhdr = spm_vol(donor);
dimg = spm_read_vols(dhdr);
hdr = spm_vol(recipient);
img = spm_read_vols(hdr);
size(img)
size(dimg)
if (size(img) ~= size(dimg)), error('Image dimensions do not match'); end;

[pth nm ext] = spm_fileparts(recipient);
hdr = dhdr;
hdr.fname = fullfile(pth, ['z' nm ext]);  
spm_write_vol(hdr,img);

