function bmp_unsharpmask(Filename, lambda, ShowFigures);
%Blurs image, uses difference between blur and orignal for edges, creates image with enhanced edges 
%  Filename: name of bitmap [optional]
%  BlurPixels: images will be smoothed with based on this FWHM
%  ShowFigures: if FALSE results saved to disk, else displayed
%Example
% bmp_unsharpmask('photo.png');
% bmp_unsharpmask('cat.jpg',0.75,-0.5,false,'N95'); %non-linear, auto-bias
% bmp_unsharpmask('cat.jpg',0.75,-0.5,true,'L95'); %linear, auto-bias
if (nargin < 1)  
   [files,pth] = uigetfile({'*.bmp;*.jpg;*.png;*.tiff;';'*.*'},'Select the Image[s]', 'MultiSelect', 'on'); 
   files = cellstr(files); %make cellstr regardless of whether user selects single or multiple images
else
    [pth,nam, ext] = fileparts(Filename);
    files = cellstr([nam, ext]); 
end;
if (nargin < 2)  
    lambda = 1.0;
end;
if (nargin < 3)
    ShowFigures = true;
end;

for i=1:size(files,2) %apply to image(s)
    nam = strvcat(deblank(files(:,i)));
    Inname = fullfile(pth, [nam]);
    Im = Imread_mat2gray_sub(Inname);
    ImOut = gradientfilt(Im);
    ImOutS = unsharp(Im,ImOut,-lambda);
    ImOutU = unsharp(Im,ImOut,lambda);
    %size(ImOut)
    if ~ShowFigures     
        imwrite(RGB_Grayscale_sub(Im),fullfile(pth, ['Original_' nam ]));
        imwrite(RGB_Grayscale_sub(ImOut+0.5),fullfile(pth, ['Edge_' nam ]));
        imwrite(RGB_Grayscale_sub(ImOutS),fullfile(pth, ['Smooth_' nam ]));
        imwrite(RGB_Grayscale_sub(ImOutU),fullfile(pth, ['UnsharpMask_' nam ]));
    end;
    
    if ShowFigures %display histogram
        figure;
        set(gcf,'color','w');
        subplot(2,2,1);
        
        
        %Im = RGB_Grayscale_sub(Im);
        %size( RGB_Grayscale_sub(Im))
        image(RGB_Grayscale_sub(Im));
        xlabel('Original');
        set(gca,'XTick',[],'YTick',[]);
        
        subplot(2,2,2);
        image(RGB_Grayscale_sub(ImOut+0.5));
        xlabel('Edges (High Frequencies)');
        set(gca,'XTick',[],'YTick',[]);
        
        subplot(2,2,3);
        image(RGB_Grayscale_sub(ImOutS));
        xlabel('Smooth (Low Frequencies)');
        set(gca,'XTick',[],'YTick',[]);
        
        subplot(2,2,4);
        image(RGB_Grayscale_sub(ImOutU));
        xlabel('Unsharp (Accentuate High Frequencies)');
        set(gca,'XTick',[],'YTick',[]);
    end;
end;

function imOut = unsharp(im, G, lambda)
% http://www.ee.columbia.edu/~madadam/4830/hw4/hw4-matlab.html#unsharp
% O = UNSHARP(IM,G,L)		Unsharp masking
% Performs unsharp masking as follows:
% 		O = IM + L*G
% where G is the gradient performed by GRADIENTFILT(IM)
% O is returned normalized to (0,1)
%G = gradientfilt(im);
imOut = im + (lambda * G);
%imOut = imnorm(imOut);
imOut = imclip(imOut);


function im = imclip(im)
im(im>1) = 1;% Clip >1 
im(im<0) = 0;% clip < 0


function imOut = imnorm(im)
%http://www.ee.columbia.edu/~madadam/4830/hw4/hw4-matlab.html#unsharp
% O = IMNORM(I)		Image normalization
% imwrite doesnt seem to do any normalizing, so i do it here.
% normalizes to the range (0,1).
imOut = im;
% normalize
imMin = min(imOut(:));
if (imMin < 0)
   imOut = imOut + abs(imMin);
end
imOut = imOut ./ max(imOut(:));

function Im=RGB_Grayscale_sub(I)
% This function transforms a grayscale image to RGB 
%min(I(:))
%max(I(:))
Im(:,:,1)=I; Im(:,:,2)=I; Im(:,:,3)=I;
Im(Im>1) = 1;% Clip >1 
Im(Im<0) = 0;% clip < 0
%img2RGB


function imOut = gradientfilt(im)
% B = GRADIENTFILT(A)		Gradient filter.
% The gradient uses the discrete laplacian;
% 	1/4 * [0 1 0
%			 1 0 1
%			 0 1 0]
% define discrete Laplacian for gradient operation. 
L = (1/4).*	[0 1 0;
   			 1 0 1;
             0 1 0];
[R,C] = size(im);
N = 3;
% pad near the edges with the symmetric extension
if mod(N,2)
   padTL = (N-1)/2;	% left and top
   padBR = (N-1)/2;	% bottom and right
else
   padTL = N/2 - 1;
   padBR = N/2;
end
offset = padTL;
paddedIm = zeros(R+N-1, C+N-1);
[padR,padC] = size(paddedIm);
% copy the main image
paddedImg((offset+1):(offset+R), (offset+1):(offset+C)) = im;
% now copy up into the symmetric extension - top, left, bottom, right
% what a pain - is there a simpler way?
paddedImg(1:padTL,(offset+1):(offset+C)) = flipud(im(1:padTL,:));
paddedImg((offset+1):(offset+R),1:padTL) = fliplr(im(:,1:padTL));
paddedImg((padR-padBR+1):padR,(offset+1):(offset+C)) = flipud(im((R-padBR+1):R,:));
paddedImg((offset+1):(offset+R),(padC-padBR+1):padC) = fliplr(im(:,(C-padBR+1):C));
% now the corners
paddedImg(1:padTL,1:padTL) = fliplr(flipud(im(1:padTL,1:padTL)));
paddedImg(1:padTL,(padC-padBR+1):padC) = fliplr(flipud(im(1:padTL,(C-padBR+1):C)));
paddedImg((padR-padBR+1):padR,(padC-padBR+1):padC) = fliplr(flipud(im((R-padBR+1):R,(C-padBR+1):C)));
paddedImg((padR-padBR+1):padR,1:padTL) = fliplr(flipud(im((R-padBR+1):R,1:padTL)));
%paddedImg	% debug
% Now apply the gradient
%imOut = zeros(padR,padC);
imOut = zeros(R,C);
for i = (offset+1):(offset+R)
   for j = (offset+1):(offset+C)
      imOut(i-offset,j-offset) = paddedImg(i,j) - sum(sum(L .* paddedImg((i-padTL):(i+padBR),(j-padTL):(j+padBR))));
   end
end

function [result] = Imread_mat2gray_sub(filename);
%this subfunction simulates
%  result =  mat2gray(double(imread(Inname)));
%without requiring the image processing toolbox
    Im =  (double(imread(filename)));
    ImSize = size(Im);
    if length(ImSize) == 2
        n = 1; %only one layer - e.g. grayscale image
    else
        n = ImSize(3); %multiple layers, e.g. color image
    end;
    for layer = 1:n
        %the following code replicates Matlab 
        % fit range so min..max is scaled to min/max..1
        scale = 1/max(max(Im(:,:,layer)));
        result(:,:,layer)= Im(:,:,layer).*scale;
        %next lines normalize values
        % so min..max are scaled to 0..1
        %mn = min(min(Im(:,:,layer))); %minumum
        %scale = 1/max(max(Im(:,:,layer))) - mn; % 1/range
        %result(:,:,layer)= (Im(:,:,layer)-mn).*scale;
    end;
    if n > 1 
        [result,U,V]=YUV_RGB_sub(result);
    end
%end Imread_mat2gray_sub

function [Y,U,V]=YUV_RGB_sub(Im)
    % This program transform RGB layers to YUV layers....
    %  By  Mohammed Mustafa Siddeq
    %  Date 25/7/2010
    Im=double(Im);
    R=Im(:,:,1); G=Im(:,:,2); B=Im(:,:,3);
    % transfom layers to YUV
    Y=((R+2*G+B)/4);
    U=R-G;
    V=B-G;
% end YUV_RGB_sub


