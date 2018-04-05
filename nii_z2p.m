function nii_z2p(fnms)
%Convert z-scores to p-values
% fnms : file name[s] of image[s] (optional)
%Examples
% nii_z2p; %use GUI
% nii_z2p('img.nii');

if ~exist('fnms','var')
	fnms = spm_select(inf,'image','Select image[s] for NaN removal'); 
end
for i=1:size(fnms,1)
    fnm = deblank(fnms(i,:));
    hdr = spm_vol(fnm);
    img = spm_read_vols(hdr);
    if size(img,4) > 1
        fprintf('%s designed for 3D images with only a single volume\n',mfilename);
    else
        [pth nm ext] = spm_fileparts(fnm);
        img(img(:) == 0) = NaN; %0 to nan
        img = spm_Ncdf(img); %convert p to z
        hdr.fname = fullfile(pth, ['p' nm ext]);  
        spm_write_vol(hdr,img);
    end
end;
