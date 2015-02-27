function bmp_yuv(Filename, Sigma, ShowFigures)
%Blurs y, u, v components of an image, saving output as bitmap with prefix
%  Filename: name of bitmap [optional]
%  Sigma: images will be smoothed with based on this FWHM
%           If Sigma is equal to zero, nearest neighbor subsampling is used
%  ShowFigures: if TRUE results displayed, else results saved
%Example
% bmp_yuv('photo.png');
% bmp_yuv('cat.jpg',0.75,-0.5,false,'N95'); %non-linear, auto-bias
% bmp_yuv('cat.jpg',0.75,-0.5,true,'L95'); %linear, auto-bias

fprintf('Demo: Humans are more sensitive to changes in brightness than hue.\n');
if (nargin < 1)  
   [files,pth] = uigetfile({'*.bmp;*.jpg;*.png;*.tiff;';'*.*'},'Select the Image[s]', 'MultiSelect', 'on'); 
   files = cellstr(files); %make cellstr regardless of whether user selects single or multiple images
else
    [pth,nam, ext] = fileparts(Filename);
    files = cellstr([nam, ext]); 
end;
if (nargin < 2)  
    Sigma = 3.0;
end;
if (nargin < 3)
    ShowFigures = true;
end;
Cutoff = ceil(3*Sigma);
for i=1:size(files,2) %apply to image(s)
    nam = char(deblank(files(:,i)));
    Inname = fullfile(pth, nam);
    Im = imread(Inname);
    ImSize = size(Im);
    if length(ImSize) == 2
        fprintf('Only able to process RGB images - not grayscale\n');
        return; %only one layer - e.g. grayscale image
    end;
    if isa(Im,'uint16')
        scale = 1/65535;
    elseif isa(Im,'uint8')
        scale = 1/255;
    else 
        fprintf('Unsupported data format\n');
        return;
    end;
    Im = double(Im) .* scale; %scale to range 0..1
    [Y,U,V]=YUV_RGB_sub(Im); %convert to Y,U,V
    if Sigma == 0
        sY = subSample(Y,2);
        sU = subSample(U,2);
        sV = subSample(V,2);
    else
        %h = fspecial('gaussian',[1,2*Cutoff+1],Sigma)% requires image pro toolbox
        h  = fspecial_sub(Cutoff, Sigma);        
        sY = conv2Sub(h,h,Y,'same'); %smoothed (blurred) Y
        sU = conv2Sub(h,h,U,'same'); %smoothed (blurred) U
        sV = conv2Sub(h,h,V,'same'); %smoothed (blurred) V
    end
    Im_sYUV=RGB_YUV_sub(sY,U,V); %image with blurred Y
    Im_YsUV=RGB_YUV_sub(Y,sU,V); %image with blurred U
    Im_YUsV=RGB_YUV_sub(Y,U,sV); %image with blurred V
    Im_YsUsV=RGB_YUV_sub(Y,sU,sV); %image with blurred U and blurred V
    if ShowFigures
        set(gcf,'color','w');

        plotImg(RGB_Grayscale_sub(Y),'Intensity (Y)',0.01,0.72);
        plotImg(RGB_GrayscaleNorm_sub(U),'R-G (U)',0.33,0.72);
        plotImg(RGB_GrayscaleNorm_sub(V),'B-G (V)',0.67,0.72);
        plotImg(RGB_Grayscale_sub(sY),'Blurred Intensity (bY)',0.01,0.39);
        plotImg(RGB_GrayscaleNorm_sub(sU),'Blurred R-G (bU)',0.33,0.39);
        plotImg(RGB_GrayscaleNorm_sub(sV),'Blurred B-G (bV)',0.67,0.39);

        plotImg(Im, 'Original', 0.01, 0.05);
        plotImg(Im_sYUV, 'Blurred Intensity (bYUV)', 0.33, 0.05);
        plotImg(Im_YsUsV, 'Blurred Color (YbUbV)', 0.67, 0.05);        
        print (gcf, '-r300', '-dpng', 'test.png'); %<- save as 300dpi , '-noui'
    else     
        imwrite(RGB_Grayscale_sub(sY),fullfile(pth, ['sY_' nam ]));
        imwrite(RGB_GrayscaleNorm_sub(sU),fullfile(pth, ['sU_' nam ]));
        imwrite(RGB_GrayscaleNorm_sub(sV),fullfile(pth, ['sV_' nam ]));
        imwrite(Im_sYUV,fullfile(pth, ['sYUV_' nam ]) );
        imwrite(Im_YsUV,fullfile(pth, ['YsUV_' nam ]) );
        imwrite(Im_YUsV,fullfile(pth, ['YUsV_' nam ]) );
        imwrite(Im_YsUsV,fullfile(pth, ['YsUsV_' nam ]) );
        Im_YUV=RGB_YUV_sub(Y,U,V); %reconstruct original
        imwrite(Im_YUV,fullfile(pth, ['YUV_' nam ]) );
    end;
end;
%end bmp_yuv()

function plotImg ( Img, Caption, X, Y)
subplot('Position',[X Y 0.3 0.28]); %width, height
image(Img);
set(gca,'XTickLabel', [],'XTick',[],'YTick',[]);
xlabel(Caption);
axis image
%end plotImg()

function B = subSample (A, N)
%Replicate every N voxels, e.g. if 3 then each 3x3 subsample has same value 
% essentially mimic JPEG sub-sampling, nearest neighbor interpolation to keep extreme values
B = A; %preallocate
yOut = 1-N;
for y = 1:size(A,1)
    if (rem(y-1,N) == 0) 
		yOut = yOut + N;
    end
    xOut = 1-N;
	for x = 1:size(A,2)
        if rem(x-1,N) == 0 
            xOut = xOut + N;
        end
		B(y,x) = A(yOut,xOut);
	end %x
end %y
%end subSample()

function H = fspecial_sub(Cutoff,Sigma)
%construct Gaussian filter
% clones Image Processing Toolbox command
%   h = fspecial('gaussian',[1,2*Cutoff+1],Sigma)
%http://www1.ynao.ac.cn/~jinhuahe/know_base/othertopics/math.htm
%   g(x) = 1/Sigma/sqrt(2*pi) * exp[-(x)^2/2/Sigma^2],
% where sigma is the standard diviation of the gaussian probability distribution. 
%   FWHM = 2*sqrt(2*log(2)) * Sigma
% in other words http://en.wikipedia.org/wiki/Full_width_at_half_maximum
%   FWHM = 2.35482004503 * Sigma
% The gaussian function can be expressed in terms of  FWHM as
%    g(x) = 2*sqrt(log(2)/pi)/FWHM * exp(-4*log(2) * (x)^2/FWHM^2)
%    h = 2*sqrt(log(2)/pi)/FWHM * exp(-4*log(2) * [-Cutoff: 1: Cutoff].^2/FWHM^2)
H = 1/Sigma/sqrt(2*pi) * exp(-(-Cutoff: 1: Cutoff).^2/2/Sigma ^2);
H = H/sum(H); %normalize results
%end fspecial_sub()

function sY = conv2Sub (h1,h2,Y,shape)
%conv2 zero pads data - here we adjust so mean of edges is zero, minimizes vignetting
if (size(Y,1) > 2) && (size(Y,2) > 2)
   mn = sum(Y(1,:)) + sum(Y(size(Y,1),:)) + sum(Y(2:size(Y,1)-1,1))+sum(Y(2:size(Y,1)-1,size(Y,2)));
   n = (2 * size(Y,1)) + (2 * size(Y,2)) - 4; 
   mn = mn / n;
else
   mn = 0; 
end
Y = Y - mn;
sY = conv2(h1,h2,Y,shape); %smoothed (blurred) Y
sY = sY + mn;
%end conv2Sub()

function Im=RGB_Grayscale_sub(I)
% This function transforms a grayscale image to RGB 
Im(:,:,1)=I; Im(:,:,2)=I; Im(:,:,3)=I;
Im(Im>1) = 1;% Clip >1 
Im(Im<0) = 0;% clip < 0
% end RGB_Grayscale_sub()

function Im=RGB_GrayscaleNorm_sub(I)
% This function transforms a grayscale image with range -1..1 to RGB with
% range 0..1
mn = min(I(:));
rng = max(I(:))-mn;
if rng ~= 0
    I = (I-mn) ./rng;
end
Im(:,:,1)=I;Im(:,:,2)=I;Im(:,:,3)=I;
%end RGB_GrayscaleX_sub()

%function Im=RGB_GrayscaleX_sub(I)
% This function transforms a grayscale image with range -1..1 to RGB with
% range 0..1
%Im(:,:,1)=(I+1)/2; Im(:,:,2)=(I+1)/2; Im(:,:,3)=(I+1)/2; 
%end RGB_GrayscaleX_sub()

function [Y,U,V]=YUV_RGB_sub(Im)
% This function transform RGB layers to YUV layers....
%  By  Mohammed Mustafa Siddeq, Date 25/7/2010
%Im=double(Im);
R=Im(:,:,1); G=Im(:,:,2); B=Im(:,:,3);
% transfom layers to YUV
Y=((R+2*G+B)/4);
U=R-G;
V=B-G;
% end YUV_RGB_sub()

function Im=RGB_YUV_sub(Y,U,V)
% This program transform YUV layers to RGB Layers in one matrix 'Im'....
%  By   Mohammed Mustafa Siddeq, Date 25/7/2010
G=((Y-(U+V)/4));
R=U+G;
B=V+G;
%next clip intensities to range 0..1
R(R<0) = 0;
G(G<0) = 0;
B(B<0) = 0;
R(R>1) = 1;
G(G>1) = 1;
B(B>1) = 1;
Im(:,:,1)=R; Im(:,:,2)=G; Im(:,:,3)=B; 
% end RGB_YUV_sub

%with regards to RGB_YUV and YUV_RGB:
% Copyright (c) 2011, Mohammed Siddeq
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.