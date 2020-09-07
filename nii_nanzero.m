function nii_nanzero(fnms)
%Voxels with zeros are replaced with NaN's, output has 'z' prefix
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
        numNan = sum(img(:)== 0);
        %img(img > 800) = 800; %example clipping extreme values
        fprintf('%s has %d voxels with zero value converted to not-a-number (NaN)\n', nm, numNan);
        if numNan > 0
            hdr.fname = fullfile(pth, ['z' nm ext]);  
            img(img == 0) = NaN;%max(img(:)); % use ~isfinite instead of isnan to replace +/-inf with zero
            %img = - img;
            spm_write_vol(hdr,img);
        end
    end
end;
