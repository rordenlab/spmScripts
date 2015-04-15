function nii_clip2bb(vols, bb, clipZ, overwrite, evenXY)
%Restrict size of image to bounding box - use after coregister
%  vols  : image(s) to reslice
%  bb    : bounding-box to preserve
%  clipZ : if true, all dimensions clipped, if false only X,Y
%          this is important for fMRI where slice timing depends on slices
%  overwrite : if true new image replaces old
%  evenXY : FSL's topup requires images have even number of rows/columns/slices,
%           If true, output image will have even number of rows/columns
% Examples
%  nii_clip2bb('T1_P001.nii',[],true); %T1 scan
%  nii_clip2bb('ASL_LM1022.nii',[],false); %4D image ASL or fMRI
%  nii_clip2bb({'APDTI_LM1021.nii','PADTI_LM1021.nii'}); %compute with AP, apply to both

if ~exist('vols','var') || isempty(vols) %no files specified
	vols = spm_select(inf,'image','Reset origin for selected image(s) (estimated from 1st)');
end
if ischar(vols), vols = cellstr(vols); end
if ~exist('bb','var') || isempty(bb) %bounding box not specified
    bb = [-78 -112 -70; 78 76 85];
end
if ~exist('clipZ','var') || isempty(clipZ) %x-clipping not specified
    clipZ = false;
end
if ~exist('evenXY','var') || isempty(evenXY) %even rows/columns specified
    evenXY = false;
end
if ~exist('overwrite','var') || isempty(overwrite) %bounding box not specified
    overwrite = false;
end
[pth,nam,ext, ~] = spm_fileparts(deblank(vols{1}));
fname = fullfile(pth,[nam ext]); %strip volume label
for i = 1: 8
	mmi = bb(1,:);
	if (mod(i,2) == 1), mmi(1) = bb(2,1); end
	if (mod(i-1,4) < 2), mmi(2) = bb(2,2); end;
	if (mod(i-1,8) < 4), mmi(3) = bb(2,3); end;
	mm(i,:) = mmi; %#ok<AGROW>
end
hdr = spm_vol([fname,',1']); %read header of 1st volume
v2m = hdr.mat; %voxel2mm transform
m2v=inv(v2m); %mm2voxel transform
for i=1:size(mm,1)
    vox(i,1:3)=mm(i,:)*m2v(1:3,1:3)' + m2v(1:3,4)'; %#ok<AGROW>
end 
mn = floor(min(vox));
mx = ceil(max(vox));
if ~clipZ %do not clip in Z-dimension - preserve number of slices
   mn(3) = 1;
   mx(3) = hdr.dim(3);
end

mxD = min(hdr.dim, mx);
mnD = max([1 1 1], mn);
if evenXY %force even number of rows and columns
	for i = 1: 2
		if mod(mxD(i)-mnD(i)+1,2) && mnD(i) > 1
			mnD(i) = mnD(i) - 1;
		end
		if mod(mxD(i)-mnD(i)+1,2) && mxD(i) < hdr.dim(i)
			mxD(i) = mxD(i) + 1;
		end
	end
	if max(mod(mxD(:)-mnD(:)+1,2))
			fprintf('Warning: odd number of rows/columns/slices may confuse topup: %s\n', fname);
	end
end %if evenXY
if (max(mnD) == 1) &&  (min(hdr.dim == mxD) == 1)
	fprintf('%s :No need to crop image\n', mfilename);
	return;
end;
%write data
vx = (mxD(1)-mnD(1)+1)*(mxD(2)-mnD(2)+1)*(mxD(3)-mnD(3)+1);
if vx <= 1, error('image not coregistered'); end;
pct = 100* vx/(hdr.dim(1)*hdr.dim(2)*hdr.dim(3));
h.dim = hdr.dim;
fprintf('%s cropping image from %dx%dx%d -> %dx%dx%d (%g%%)\n', mfilename, h.dim(1), h.dim(2), h.dim(3), mxD(1)-mnD(1)+1,mxD(2)-mnD(2)+1, mxD(3)-mnD(3)+1, pct);
for v = 1 : numel(vols) %apply parameters from first session to others
    [pth,nam,ext, ~] = spm_fileparts(deblank(vols{v}));
    fname = fullfile(pth,[nam ext]); %strip volume label
    hdr = spm_vol([fname,',1']); %read header of 1st volume
    if (h.dim(1) ~= hdr.dim(1)) || (h.dim(2) ~= hdr.dim(2)) || (h.dim(3) ~= hdr.dim(3))
        error('%s error: Image dimensions do not match %dx%dx%d ~= %dx%dx%d %s %s', mfilename, ...
            h.dim(1),h.dim(2),h.dim(3), hdr.dim(1),hdr.dim(2),hdr.dim(3), deblank(vols{1}), deblank(vols{v}));
        
    end
end
for v = 1 : numel(vols) %apply parameters from first session to others
    [pth,nam,ext, ~] = spm_fileparts(deblank(vols{v}));
    fname = fullfile(pth,[nam ext]); %strip volume label
    hdr = spm_vol(fname); %read header - this time 4D if specified
    img = spm_read_vols(hdr); %load image
    hdr = spm_vol([fname,',1']); %read header of 1st volume
    if (h.dim(1) ~= hdr.dim(1)) || (h.dim(2) ~= hdr.dim(2)) || (h.dim(3) ~= hdr.dim(3))
        error('%s error: Image dimensions do not match %s %s', mfilename, deblank(vols{1}), deblank(vols{v}));
        
    end
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
end
%spm_write_vol(hdr,img);