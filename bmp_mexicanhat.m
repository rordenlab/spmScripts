function bmp_mexicanhat(Filename,  ShowFigures)
%Blurs image, uses difference between blur and orignal for edges, creates image with enhanced edges 
%  Filename: name of bitmap [optional]
%  ShowFigures: if FALSE results saved to disk, else displayed
%Example
% bmp_unsharpmask('photo.png');
% bmp_unsharpmask('cat.jpg',0.75,-0.5,false,'N95'); %non-linear, auto-bias
% bmp_unsharpmask('cat.jpg',0.75,-0.5,true,'L95'); %linear, auto-bias
if ~exist('Filename','var')  
   [files,pth] = uigetfile({'*.bmp;*.jpg;*.png;*.tiff;';'*.*'},'Select the Image[s]', 'MultiSelect', 'on'); 
   files = cellstr(files); %make cellstr regardless of whether user selects single or multiple images
else
    [pth,nam, ext] = fileparts(Filename);
    files = cellstr([nam, ext]); 
end;
if ~exist('ShowFigures','var')
    ShowFigures = true;
end;
for i=1:size(files,2) %apply to image(s)
    nam = strvcat(deblank(files(:,i))); %#ok<REMFF1>
    Inname = fullfile(pth, nam);
    Im = Imread_mat2gray_sub(Inname);
    ImOut = logImgSub(Im, 32,16);
    ImOut2 = logImgSub(Im, 32,0.1);
    ImOut = normSub(ImOut);
    ImOut2 = normSub(ImOut2);
    %size(ImOut)
    if ~ShowFigures     
        imwrite(RGB_Grayscale_sub(Im),fullfile(pth, ['Original_' nam ]));
        imwrite(RGB_Grayscale_sub(ImOut),fullfile(pth, ['Edge_' nam ]));
        imwrite(RGB_Grayscale_sub(ImOut2),fullfile(pth, ['Edge2_' nam ]));
    end;
    if ShowFigures %display histogram
        figure;
        set(gcf,'color','w');
        subplot(1,3,1);
        image(RGB_Grayscale_sub(Im));
        xlabel('Original');
        set(gca,'XTick',[],'YTick',[]); 
        subplot(1,3,2);
        image(RGB_Grayscale_sub(ImOut));
        xlabel('Edges (Medium Frequencies)');
        set(gca,'XTick',[],'YTick',[]);
        subplot(1,3,3);
        image(RGB_Grayscale_sub(ImOut2));
        xlabel('Edges (High Frequencies)');
        set(gca,'XTick',[],'YTick',[]);
    end;
end;

function imOut = normSub(im)
% normalizes to the range (0,1).
imOut = im;
imOut = imOut - min(imOut(:));
if max(imOut(:)) == 0
    return
end
imOut = imOut ./ max(imOut(:));
%imOut = sqrt(imOut);
%imOut = imOut.^2;
%end imnorm

function imOut = logImgSub(im, hsize, sigma)
%apply Laplacian of Gaussian (LoG) to image
%mexhat = [0 1 0; 1 -4 1; 0 1 0];
mexhat = logSub(hsize, sigma);
imOut = conv2(im,mexhat,'same'); 
%end gradientfilt()

function f = logSub(hsize,sigma)
%estimate Laplacian of Gaussian (LoG) 
%http://nesl.ee.ucla.edu/fw/chenni/Chenni_desktop/Documents/gumstix_ubuntu/2.0%20GB%20Filesystem/usr/share/octave/packages/3.0/image-1.0.8/fspecial.m
if ~exist('hsize','var')
    hsize = [5, 5];
else
    if (length(hsize(:)) == 1)
      hsize = [hsize, hsize];
    end
end
% Get sigma
if ~exist('sigma', 'var')
    sigma = 0.5;
end
%Compute the filter
h1 = hsize(1)-1; h2 = hsize(2)-1; 
[x y] = meshgrid(0:h2, 0:h1);
x = x-h2/2; 
y = y-h1/2;
gauss = exp( -( x.^2 + y.^2 ) / (2*sigma^2) );
f = ( (x.^2 + y.^2 - 2*sigma^2).*gauss )/( 2*pi*sigma^6*sum(gauss(:)) );
f = f - mean(f(:));
%end logSub()

function Im=RGB_Grayscale_sub(I)
% This function transforms a grayscale image to RGB 
%min(I(:))
%max(I(:))
Im(:,:,1)=I; Im(:,:,2)=I; Im(:,:,3)=I;
Im(Im>1) = 1;% Clip >1 
Im(Im<0) = 0;% clip < 0
%end RGB_Grayscale_sub()

function [result] = Imread_mat2gray_sub(filename)
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
    result = Im; %pre-allocate output
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
        [result,~,~]=YUV_RGB_sub(result);
    end
%end Imread_mat2gray_sub()

function [Y,U,V]=YUV_RGB_sub(Im)
% This function transform RGB layers to YUV layers....
%  By  Mohammed Mustafa Siddeq, Date 25/7/2010
%Im=double(Im);
R=Im(:,:,1); G=Im(:,:,2); B=Im(:,:,3);
% transfom layers to YUV
Y=((R+2*G+B)/4);
U=R-G;
V=B-G;
%end YUV_RGB_sub()
        
        
        

