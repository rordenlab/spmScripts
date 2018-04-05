function nii_p2z(fnms)
%Convert p values to z scores
% fnms : file name[s] of image[s] (optional)
%Examples
% nii_p2z; %use GUI
% nii_p2z('img.nii');

if ~exist('fnms','var')
	fnms = spm_select(inf,'image','Select image[s] for NaN removal'); 
end
for i=1:size(fnms,1)
    fnm = deblank(fnms(i,:));
    hdr = spm_vol(fnm);
    img = spm_read_vols(hdr);
    if (max(img(:)) > 1) || (min(img(:)) < 0) 
        error('p-values should range from 0..1 not %g..%g', min(img(:)), max(img(:)) );
    end;
    if size(img,4) > 1
        fprintf('%s designed for 3D images with only a single volume\n',mfilename);
    else
        [pth nm ext] = spm_fileparts(fnm);
        img(img(:) == 0) = NaN; %0 to nan
        img = spm_invNcdf(img); %convert p to z
        %img = -img; %reverse polarity: spm_invNcdf([0.05,0.95]) yields [-1.6449, 1.6449]
        hdr.fname = fullfile(pth, ['z' nm ext]);  
        spm_write_vol(hdr,img);
    end
end;
