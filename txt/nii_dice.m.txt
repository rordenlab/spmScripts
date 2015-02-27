function d = nii_dice(im1, im2, binarize, im1inten, im2inten);
%computes Dice's coefficient (see also nii_jaccard)
%http://en.wikipedia.org/wiki/Dice's_coefficient
%http://sve.loni.ucla.edu/instructions/metrics/dice/
% input
%  im1: image
%  im2: image
%  binarize [optional]: if true voxels are treted as 1 or 0
% output
%  d: Dice's coefficient
%d =nii_dice('c1head_1.nii','mask_gray.nii')

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


inter = sum(sum(sum(min(i1,i2)) ) ); 
s1 = sum(i1(:));
s2 = sum(i2(:));

%in terms of intersection between images x and y (iX_Y) versus sum of X (sX) and sum of Y (sY), d = (2*iX_Y)/(s1+s2)
%in terms of false positive(fp), false negative(fn), true negative(tn), and true positive (tp) counts, d = 2*TP/((FP+TP)+(TP +FN))
d = (2*inter) / (s1 + s2);
fprintf('Dice coefficient of %s and %s is %f\n',im1.fname, im2.fname, d);

