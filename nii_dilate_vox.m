function img = nii_dilate_vox(fnms, dilatevox)
%dilate binary volume
% dilatevox: Number of voxels to grow or shrink object
%   if dilatevox = 0: no change to object
%   if dilatevox = 2: output will be dilated 2 voxels larger than input
%By default, images are saved to disk
%  However, if img is used as an output no data is saved to disk
%  if output is used, fnms should only specify one image
%Example:
%  nii_dilate_vox('32mmBall.nii',2); %dilate

if ~exist('fnms','var') %no files
    fnms = spm_select(inf,'image','Select images to dilate');
end
if ~exist('dilatevox','var') %no FWHM specified
    dilatevox = inputdlg({'Dilation in voxels'},'Enter integer', [1],{'2'});
    dilatevox = str2double(dilatevox{1});
end
if dilatevox == 0
    error('dilatevox should be a non-zero integer');
end
for i=1:size(fnms,1)
  fnm = deblank(fnms(i,:));
  [pth,nam,ext] = spm_fileparts(fnm);
  src = fullfile(pth,[nam ext]);
  hdr = spm_vol(src);
  img = spm_read_vols(hdr);
  img(isnan(img(:))) = 0; %remove NAN
  img(img(:) ~= 0) = 1; %binarize
  if dilatevox < 0
    img = 1.0 - img;
    df = distance_field3D (img);
    %hdr.fname = fullfile(pth,['df',  nam, ext]);
    %hdr.dt(1)=16;%save binary data as bytes: uint8=2; int16=4; int32=8; float32=16; float64=64
    %hdr.pinfo(1:2) = [1 0]; %scl_slope, scl_inter
    %spm_write_vol(hdr, df);
    %return;
    img = 1.0 - (df <= abs(dilatevox));
  else
    df = distance_field3D (img);
    img = (df <= dilatevox)+0;
  end
  if nargout > 0 
    return;
  end
  hdr.fname = fullfile(pth,['d',  nam, ext]);
  hdr.dt(1)=2;%save binary data as bytes: uint8=2; int16=4; int32=8; float32=16; float64=64
  hdr.pinfo(1:2) = [1 0]; %scl_slope, scl_inter
  spm_write_vol(hdr, img);  
end
%end nii_dilate_vox()

function df = distance_field3D (img)
% compute distance field for a three-dimensional binary image.
% each voxel in the distance field is assigned to the distance to the
% nearest "true" voxel of the binary image.
% the Python code for two-dimensional image is at Philip Rideout's blog:
% https://prideout.net/blog/distance_fields/
% the pseudocode for one-dimensional vector is here:
% Felzenszwalb, P. F., & Huttenlocher, D. P. (2012). Distance transforms of sampled functions. Theory of computing, 8(1), 415-428.
% Matlab translation by Grigori Yourganov
img (img == 0) = Inf;
img (img == 1) = 0;
for i = 1:size (img, 2)
    for j = 1:size (img, 3)
        img1 (:, i, j) = horizontal_pass (img (:, i, j));
    end
end
img2 = permute (img1, [2 1 3]);
for i = 1:size (img2, 2)
    for j = 1:size (img2, 3)
        img3 (:, i, j) = horizontal_pass (img2 (:, i, j));
    end
end
img2 = permute (img3, [3 2 1]);
clear img3
for i = 1:size (img2, 2)
    for j = 1:size (img2, 3)
        img3 (:, i, j) = horizontal_pass (img2 (:, i, j));
    end
end
img3 = permute (img3, [2 3 1]);
df = sqrt (img3);


function df = horizontal_pass (single_row)
[hull_vertices, hull_intersections] = find_hull_parabolas (single_row);
df = march_parabolas (single_row, hull_vertices, hull_intersections);


function [hull_vertices, hull_intersections] = find_hull_parabolas (single_row)
k = 1;
v(1) = 1;
z(1) = -Inf;
z(2) = Inf;
for q = 2:length (single_row)
    s = intersect_parabolas (single_row, q, v(k));
    while s <= z(k) && k > 1
        k = k - 1;
        s = intersect_parabolas (single_row, q, v(k));
    end
    k = k + 1;
    v(k) = q;
    z(k) = s;
    z(k+1) = Inf;
end
hull_vertices = v;
hull_intersections = z;

function s = intersect_parabolas (single_row, p, q)
s = ((single_row (p) + p*p) - (single_row (q) + q*q)) / (2*p - 2*q);
if isnan (s)
    s = Inf;
end

function df = march_parabolas (single_row, hull_vertices, hull_intersections)
d = single_row;
v = hull_vertices;
z = hull_intersections;
k = 1;
for q = 1:length (d)
    while z(k+1) < q
        k = k + 1;
    end
    dx = q - v(k);
    single_row(q) = dx*dx + single_row(v(k));
end
df = single_row;
