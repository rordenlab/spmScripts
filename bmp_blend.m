function bmp_blend (A,B,Frac,Outname)
%Creates weighted average of image A and B
% Frac determines weight of image A: 0.75 yields 75%A 25%B 
%  A,B: names of bitmaps to be blended [optional]
%  Frac: Weight for image A: a value 0..1 [optional]
%  Outname: name for blended image [optional]
%Example
% bmp_blend('cat.png','dog.png',0.25);

if (nargin < 1)  %A not specified
   [A,Apth] = uigetfile({'*.bmp;*.jpg;*.png;*.tiff;';'*.*'},'Select first image');
   A = [Apth, A];
end;
if (nargin < 2)  %B not specified
   [B,Bpth] = uigetfile({'*.bmp;*.jpg;*.png;*.tiff;';'*.*'},'Select second image');
   B = [Bpth, B];
end;
if (nargin < 3) 
    Frac = 0.5; 
end; %Frac not specified
if (nargin < 4) %outname not specified
    [pth,nam, ext] = fileparts(A);
	Outname = fullfile(pth, ['b' nam ext]);
end;
if ((Frac > 1) || (Frac < 0)) 
    Frac = 0.5; 
end; %Frac must be in range 0..1

%process images
Resid = 1 - Frac;
Ia =  imread(A);
Ib =  imread(B);
assert(all(size(Ia)==size(Ib)),'bmp_blend error: images have different sizes');
Iblend(:,:,:) = Frac*Ia(:,:,:) + Resid*Ib(:,:,:); 
imwrite(Iblend,Outname);
%imshow(Iblend) % <- requires image processing toolbox

