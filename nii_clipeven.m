function nii_clipeven(vols, overwrite)
%Make sure slices have even number of rows/columns (required by topup)
%  vols  : image(s) to reslice
% Examples
%  nii_clipeven('T1_P001.nii'); 

if ~exist('vols','var') || isempty(vols) %no files specified
	vols = spm_select(inf,'image','Select images to reslice');
end
if ~exist('overwrite','var') || isempty(overwrite) %bounding box not specified
    overwrite = false;
end
if ischar(vols), vols = cellstr(vols); end

[pth,nam,ext, ~] = spm_fileparts(deblank(vols{1}));
fname = fullfile(pth,[nam ext]); %strip volume label
hdr = spm_vol([fname,',1']); %read header of 1st volume
if ~mod(hdr.dim(1),2) && ~mod(hdr.dim(2),2) 
    fprintf('%s :No need to crop image (even numer of rows and columns)\n', mfilename);
    return;
end;
for i = 1:3
    mxD(i) = hdr.dim(i);
    mnD(i) = 1;
end
if mod(hdr.dim(1),2)
    mxD(1) = mxD(1) - 1;
end
if mod(hdr.dim(2),2)
    mxD(2) = mxD(2) - 1;
end
vx = (mxD(1)-mnD(1)+1)*(mxD(2)-mnD(2)+1)*(mxD(3)-mnD(3)+1);
if vx <= 1, error('Can not crop volumes with only one row/column'); end;
pct = 100* vx/(hdr.dim(1)*hdr.dim(2)*hdr.dim(3));
fprintf('%s cropping image from %dx%dx%d -> %dx%dx%d (%g%%)\n', mfilename, hdr.dim(1), hdr.dim(2), hdr.dim(3), mxD(1)-mnD(1)+1,mxD(2)-mnD(2)+1, mxD(3)-mnD(3)+1, pct);
v2m = hdr.mat; %voxel2mm transform
m2v=inv(v2m); %mm2voxel transform
for v = 1 : numel(vols) %apply parameters from first session to others
    [pth,nam,ext, ~] = spm_fileparts(deblank(vols{v}));
    fname = fullfile(pth,[nam ext]); %strip volume label
    hdr = spm_vol([fname]); %read header - this time 4D if specified
    img = spm_read_vols(hdr); %load image
    hdr = spm_vol([fname,',1']); %read header of 1st volume
    img = img(mnD(1):mxD(1), mnD(2):mxD(2), mnD(3):mxD(3), :); %clip image dimensions
    origin= mnD*v2m(1:3,1:3)' + v2m(1:3,4)';
    hdr.mat(1:3,4) = origin;
    hdr.dim = size(img); 
    hdr.dim = hdr.dim(1:3); %for 4D volumes, treat each volume separately
    if ~overwrite
        hdr.fname = fullfile(pth,['c' nam ext]);
    else
        delete(fname);
    end
    for vol=1:size(img,4)
        hdr.n(1)=vol;
        spm_write_vol(hdr,img(:, :, :, vol));
    end;
end %each volumne
%end nii_clipeven