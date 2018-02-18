function [outhdr, outimg] = nii_reslice_target(inhdr, inimg, tarhdr, interp) 
%Reslice input image to match dimensions of target image (either to disk or memory)
%  inhdr: image to reslice- either filename of NIfTI header or loaded NIfTI header structure
%  inimg: (optonal) NIfTI image data (only if inhdr is a structure)
%  tarhdr: image to match- either filename of NIfTI header or loaded NIfTI header structure
%  interp: (optional) if 0 nearest neighbor interpolation, else trilinear interpolation (default)
% Outputs: if not specified, resliced image saved to disk, else returns resliced header and image
%Chris Rorden (2014) see John Ashburner and Ged Ridgway's reorient.m 
% http://opensource.org/licenses/BSD-2-Clause
%Examples
%   nii_reslice_target_thresh('spmT_5.nii', '','mni152_2009_256.nii')
% Next example shows how to use this without writing to disk
%  inhdr = spm_vol('T1_LM1001.nii'); %load header
%  inimg = spm_read_vols(inhdr); %load volume
%  tarhdr = spm_vol('jhu.nii'); %load header
%  [outhdr,outimg] = nii_reslice_target_thresh(inhdr, inimg, tarhdr); %resize in memory
%  spm_write_vol(outhdr,outimg); %save resized image
% Next example: specify target dimensions without explicitly providing a full header
%  tarhdr.mat = [-1 0 0 79; 0 1 0 -113; 0 0 1 -51; 0 0 0 1];
%  tarhdr.dim = [157 189 136];
%  nii_reslice_target_thresh('qCBV.nii','',tarhdr);

if ~exist('inhdr','var')
    inhdr = spm_select(1,'image','Select source image that will be resliced');
end
if ~exist('tarhdr','var')
    tarhdr = spm_select(1,'image','Select target image (source will be resliced to match target)');
end
if ~exist('interp','var')
    interp = 1;%linear
end
if ~isnumeric(interp), interp = interp + 0; end; %change false/true to 0/1
if ~isstruct(tarhdr)
    tarhdr = spm_vol(tarhdr); %load target header
end
imgdim = tarhdr.dim(1:3);
if ~isstruct(inhdr)
    inhdr = spm_vol(inhdr); %load input header
    inhdr = inhdr(1); %if 4D, only process first volume
	inimg = spm_read_vols(inhdr); %load input image
end
if size(inimg,4) > 1
    inhdr = inhdr(1);
    inimg = inimg(:,:,:,1);
    fprintf('%s warning: only reslicing first volume of 4D image %s\n', mfilename, inhdr.fname);
end
mask = (inimg == 0) + (isnan(inimg));
isMasked = sum(mask(:)) > 0;
if isMasked
    fprintf('Assuming %d voxels are masked %d=zero, %d=NaN\n', sum(mask(:)), sum(inimg(:) == 0), sum(isnan(inimg(:))) );
    if (interp == 0)
        warning('Masking has no influence on nearest neighbor interpolation');
    end
else
    warning('Image does not appear thresholded');
end;
negThresh = max(inimg(inimg(:) < 0));
posThresh = min(inimg(inimg(:) > 0));
if ~isempty(negThresh) && ~isempty(posThresh) 
   error('Not simple threshold: image has both posiitve and negative values.'); 
end
if isempty(negThresh)
    thresh = posThresh;
else
   thresh = negThresh; 
end
fprintf('Threshold %g\n',thresh);
if interp > 1
   error('Not designed for higher-order interpolation'); 
end
mask = 1 - mask;
inimg(isnan(inimg)) = 0; % convert nan to zero - crucial or major contamination
outhdr            = inhdr;
[pth,nam,ext] = fileparts(outhdr.fname);
outhdr.fname      = fullfile(pth,['x' nam ext]);
outhdr.dim(1:3)   = imgdim(1:3);
outhdr.mat        = tarhdr.mat;
if isequal(inhdr.mat, outhdr.mat) && isequal(outhdr.dim(1:3), inhdr.dim(1:3))
    fprintf('%s no need to reslice: input image already aligned to target image\n', mfilename);
    outimg = inimg;
else %if reslicing is required
    outimg = zeros(outhdr.dim(1:3));
    for i = 1:imgdim(3)
        M = inv(spm_matrix([0 0 -i])*inv(outhdr.mat)*inhdr.mat);
        outslice = spm_slice_vol(inimg, M, imgdim(1:2), interp); % (1=linear interp)
        maskslice = spm_slice_vol(mask, M, imgdim(1:2), interp); % (1=linear interp)
        %if false
            maskslice = double(maskslice >= 0.5);
            outslice(maskslice == 0) = nan;
            maskslice = maskslice * thresh;
            if thresh > 0
                outslice = max(maskslice, outslice);
            else
                outslice = min(maskslice, outslice);
            end
        %else
        %    maskslice(maskslice < 0.5) = NaN; %without this, things get dilated
        %    outslice = outslice ./ maskslice;
        %end
        outimg(:,:,i) = outslice;
    end
end
if nargout < 2
    outhdr = spm_create_vol(outhdr); %save header to disk
    spm_write_vol(outhdr,outimg); %save image to disk
end
%end nii_reslice_target()