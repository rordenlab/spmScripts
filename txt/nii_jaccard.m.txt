function j = nii_jaccard(im1, im2, binarize, im1inten, im2inten);
%computes Jaccard Similarity Index
%http://en.wikipedia.org/wiki/Jaccard_index
% input
%  im1: image
%  im2: image
%  binarize [optional]: if true Tanimoto's Similarity is computed 
% output
%  j: Jaccard Index (intersection vs union)
%j =nii_jaccard('c1head_1.nii','mask_gray.nii')

if nargin <1 %no im1
 im1 = spm_select(1,'image','Select 1st image');
end;
if nargin <2 %no im2
 im2 = spm_select(1,'image','Select 2nd image');
end;
if nargin <3 %no binarize
 binarize = false;
end;

if ischar(im1), im1 = spm_vol(im1); end;
if ischar(im2), im2 = spm_vol(im2); end;

i1 = spm_read_vols(im1);
i2 = spm_read_vols(im2);

if nargin > 3 %filter im1
    i1 = i1 == im1inten;
    i1 = uint8(i1);
    fprintf('A total of %d voxels in %s have an intensity of %f\n',sum(i1(:)), im1.fname, im1inten);
end;
	
if nargin > 4 %filter im2
    i2 = i2 == im2inten;
    i2 = uint8(i2);
    fprintf('A total of %d voxels in %s have an intensity of %f\n',sum(i2(:)), im2.fname, im2inten);
end;

if binarize  %scale intensity to 0..1, values less than 0.5 set to zero, others set to 1
    mx=max(i1(:));
    mn=min(i1(:));
    i1 = (i1-mn)/mx;
    i1 = uint8(round(i1));
    mx=max(i2(:));
    mn=min(i2(:));
    i2 = (i2-mn)/mx;
    i2 = uint8(round(i2));
end;

%we can compute intersection/union using vectors...
inter = sum(sum(sum(min(i1,i2)) ) ); %instersection: result only 1 when i1=1 and i2=1
union = sum(sum(sum(max(i1,i2)) ) );  %instersection: result only 1 when i1=1 and i2=1
%or we could compute the same values with for loops...
% inter = 0;
% union = 0;
% for i=1:size(i1,3),
%         p1       = double(i1(:,:,i));
%         p2       = double(i2(:,:,i));
%         a= max(p1,p2);
%         union = union + sum(a(:));
%         a= min(p1,p2);
%         inter = inter + sum(a(:));
% end;

j = inter/union;
fprintf('Jaccard of %s and %s is %f\n',im1.fname, im2.fname, j);

