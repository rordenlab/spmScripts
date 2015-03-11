function nii_neper(fnms)
%Convert FA image using Neper scale
% fnms : file name[s] of image[s] (optional)
%Examples
% nii_neper; %use GUI
% nii_neper('DTI_P125_FA.nii');


if ~exist('fnms','var')
	fnms = spm_select(inf,'image','Select image[s] for NaN removal'); 
end
for i=1:size(fnms,1)
    fnm = deblank(fnms(i,:));
    hdr = spm_vol(fnm);
    img = spm_read_vols(hdr);
    [pth nm ext] = spm_fileparts(fnm);
    if min(img(:)) < 0, error('Image intensity must be positive'); end;
    if min(img(:)) == max(img(:)), error('Image has no variability'); end;
    img = img - min(img(:)); %scale from zero
    img = img/max(img(:)); %scale from 0..1
    img = power(img, 0.5);
    hdr.fname = fullfile(pth, ['n' nm ext]);  
    spm_write_vol(hdr,img);
end;
