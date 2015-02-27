function nii_flipz(fnms)
%Voxels with NaN's are replaced with zero, output has 'z' prefix
% fnms : file name[s] of image[s] (optional)
%Notes
% For slower (2D) implementation see http://blogs.warwick.ac.uk/nichols/entry/zero_nans_in/
%Examples
% nii_zeronan; %use GUI
% nii_zeronan('img.nii');
% nii_zeronan(strvcat('a.nii','b.nii'));

if ~exist('fnms','var')
	fnms = spm_select(inf,'image','Select image[s] for NaN removal'); 
end
for i=1:size(fnms,1)
    fnm = deblank(fnms(i,:));
    [pth nm ext] = spm_fileparts(fnm); %remove volume number
    hdr = spm_vol(fullfile(pth, [nm, ext]) ); %without volue
    img = spm_read_vols(hdr);
    slices = size(img,3);
    if slices > 1
        hdr = hdr(1);
        hdr.fname = fullfile(pth, ['z' nm ext]); 
        imgout = img(:,:,:,1);
        for v = 1: size(img,4) % for each volume
            imgin = img(:,:,:,v);
            for z = 1: slices % for each slice (z dimension)
                imgout(:,:,z) = imgin(:,:,1+slices-z);
            end
            hdr.n = [v 1];
            spm_write_vol(hdr,imgout);            
        end


    end
end
