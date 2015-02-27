function nii_zeronan(fnms)
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
    hdr = spm_vol(fnm);
    img = spm_read_vols(hdr);
    if size(img,4) > 1
        fprintf('%s designed for 3D images with only a single volume\n',mfilename);
    else
    	[pth nm ext] = spm_fileparts(fnm);
        numNan = sum(isnan(img(:)));
        fprintf('%s has %d voxels with not-a-number (NaN) intensities\n', nm, numNan);
        if numNan > 0
            hdr.fname = fullfile(pth, ['z' nm ext]);  
            img(isnan(img)) = 0;%max(img(:)); % use ~isfinite instead of isnan to replace +/-inf with zero
            spm_write_vol(hdr,img);
        end
    end
end;
