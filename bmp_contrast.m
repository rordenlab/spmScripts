function bmp_contrast(Filename, Gain, Bias, Smooth, Linear,  Prefix, ShowFigures)
%Adjusts contrast of image, saving output as bitmap with prefix
%  Filename: name of bitmap [optional]
%  Gain: intensity amplifaction 0..1: 0.1 = low contrast, 0.5= unchanged, 0.9 high contrast
%  Bias: intensity offset 0..1: 0.1 = darker, 0.5= unchanged, 0.9  brighter
%        Special value: if bias is <0, then adaptive so image does not get darker/brighter
%  Smooth: Gaussian blur (0=none, 1=some)
%  Linear: true/false - select between linear transform or nonlinear
%  Prefix: string appended to output name, e.g. if 'e' then cat.jpg ->  ecat.jpg
%  ShowFigures: Determine whether histogram and function are graphed
%Example
% bmp_contrast('dog.png');
% bmp_contrast('cat.jpg',0.75,-0.5,0, false,'N95'); %non-linear, auto-bias
% bmp_contrast('cat.jpg',0.75,-0.5,0,true,'L95'); %linear, auto-bias

%Alternatives
% If you have the Matlab Image Processing Toolbox, you may want to consider contrast_rgb or nlfilter
%   http://people.sc.fsu.edu/~jburkardt/m_src/image_contrast/image_contrast.html
%   http://www.mathworks.com/products/image/examples.html?file=/products/demos/shipping/images/ipexcontrast.html
if (nargin < 1)  
   [files,pth] = uigetfile({'*.bmp;*.jpg;*.png;*.tiff;';'*.*'},'Select the Image[s]', 'MultiSelect', 'on'); 
   files = cellstr(files); %make cellstr regardless of whether user selects single or multiple images
else
    [pth,nam, ext] = fileparts(Filename);
    files = cellstr([nam, ext]); 
end;
if (nargin < 5)  
    %Gain = 0.7;
    prompt = {'Enter contrast (0..1; <0.5 lower, >0.5 higher):','Enter brightness (0..1; <0.5 darker, >0.5 brighter, negative for automatic):','Enter blur (0=none, 1=a little)','Linear (1) or Nonlinear (0) transform'};
    dlg_title = 'Values for adjusting the image(s)';
    num_lines = 1;
    def = {'0.8','-0.5','0','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    Gain = str2double(answer{1});
    Bias = str2double(answer{2});
    Smooth = str2double(answer{3});
    Linear = str2double(answer{4});

end;
if (nargin < 6)  
    Prefix = 'c';
end;
if (nargin < 7)
    ShowFigures = true;
end;

%apply to image(s)
for i=1:size(files,2)
    nam = strvcat(deblank(files(:,i))); %#ok<REMFF1>
    Inname = fullfile(pth, nam);
    Outname = fullfile(pth, [Prefix nam ]);
    %Im = Imread_mat2gray_sub(Inname); 
    Im = imread(Inname);
    ImSize = size(Im);
    %determine layers: grayscale =1, red/green/blue =3
    if length(ImSize) == 2
        fprintf('Only able to process RGB images - not grayscale\n');
        return; %only one layer - e.g. grayscale image
    end;
    %scale to range 0..1
    if isa(Im,'uint16')
        scale = 1/65535;
    elseif isa(Im,'uint8')
        scale = 1/255;
    else 
        fprintf('Unsupported data format\n');
        return;
    end;
    Im = double(Im) .* scale;
    %convert to Y,U,V
    [Y,U,V]=YUV_RGB_sub(Im);
    %make a histogram - much more rapid for autobalance, and creates useful graph
    wid = size(Y,1);
    ht = size(Y,2);
    Y1d = reshape(Y,ht*wid,1);
    x = 0:1/255:1; 
    h = hist(Y1d,x);
    if ShowFigures %display histogram
        figure;
        hist(Y1d,x);
        
    end; %if ShowFigures
    %compute mean
    delta=1e-6; % a very small number, so zeros do not cause problems...
    % logMeanY= logMean_sub (h); %you can compute from either histogram or based on each voxel
    logMeanY=exp(mean(mean(log(Y+delta))));%overall luminance: log average of image
    %meanY = mean(mean(Y)) ;
    fprintf(' Before transform: log Mean intensity %f\n', logMeanY);
    %compute transform
    if Bias < 0, %autobalance
        bestFit = Inf;
        hR = h;
        for ci = 0:255,
        	out = makeTransform_sub (Gain, ci/256, Linear);
            for x = 1:256, %initialize - several input intensities may be set to same output
                hR(x) = 0;
            end; 
            for x = 1:256,
                pos = round(out(x)*255) + 1; 
                hR(pos) = hR(pos) + h(x);
            end; %for x
            logMeanYr= logMean_sub (hR);
            fit = abs(logMeanYr-logMeanY);
            if (fit < bestFit) 
                Bias = ci/256;
                bestFit = fit;
            end;
        end;
        fprintf(' Auotbalance set the to bias %f\n',Bias);
    end; %if Bias < 0: autobalance 
     out = makeTransform_sub (Gain, Bias, Linear);
     if ShowFigures
        plotTransform_sub (out);
     end; %if ShowFigures
     %apply transform
     for y = 1:ImSize(2),
         for x = 1:ImSize(1),
                  Y(x,y) = out(1+round(255*Y(x,y)));
         end;
     end;
     if Smooth > 0
         preMeanY = mean(Y(:));
         %preMeanY=exp(mean(mean(log(Y+delta))));%overall luminance: log average of image
         Y = GaussianFilter(Y, 3, Smooth);
         meanY = mean(Y(:));
         %meanY=exp(mean(mean(log(Y+delta))));%overall luminance: log average of image
         Y = Y + (preMeanY - meanY);
         
     end
     %report effects
     logMeanY=exp(mean(mean(log(Y+delta))));%overall luminance: log average of image
     fprintf(' After transform: log Mean intensity %f\n', logMeanY);
     ImAdjusted = RGB_YUV_sub(Y,U,V);
     imwrite(ImAdjusted,Outname);
end;
%end bmp_contrast()

function Filtered = GaussianFilter(ImageData, hsize, sigma)
% http://stackoverflow.com/questions/13193248/how-to-make-a-gaussian-filter-in-matlab
%Get the result of Gaussian
filter_ = Gaussian2D(hsize, sigma);
Filtered = conv2(ImageData, filter_, 'same');
%%check image
%[r, c] = size(ImageData);
%Filtered = zeros(r, c);    
% for i=1:r
%     for j=1:c
%         for k=1:hsize
%             for m=1:hsize
%                     Filtered =  Filtered + ImageData(i,j).*filter_(k,m);    
%             end
%         end
%     end
% end
%end GaussianFilter()

function h =  Gaussian2D(hsize, sigma)
n1 = hsize;
n2 = hsize;
for i = 1 : n2 
    for j = 1 : n1
        % size is 10;
        % -5<center<5 area is covered.
        c = [j-(n1+1)/2 i-(n2+1)/2]';                
        % A product of both axes is 2D Gaussian filtering
        h(i,j) = Gauss(c(1), sigma)*Gauss(c(2), sigma);        
    end
end
%end Gaussian2D()

function Gaussian_filtered = Gauss(image_x, sigma)
% for single axis
% http://en.wikipedia.org/wiki/Gaussian_filter
Gaussian_filtered = exp(-image_x^2/(2*sigma^2)) / (sigma*sqrt(2*pi)); 
%end Gauss()


function [out]= logMean_sub (histo8bit)
delta=1e-6; % a very small number, so zeros do not cause problems...
%logMeanY=exp(mean(mean(log(Y+delta))));   
n = sum(histo8bit); %number of pixels
logV = 0.0;
for i = 1:256
    logV = logV + (histo8bit(i)*log( ((i-1)/255)+delta)); 
end;
out = exp(logV/n);


function plotTransform_sub (out)
%OPTIONAL: plot lookup table
in = 0:1/255:1;
figure; %to save this image rather than overwrite
p =plot(in, out);
axis([0 1 0 1]);
xlabel('Input Intensity');
ylabel('Output intensity');
%legend('Gray Matter');
title( 'Contrast Correction');
set(p,'LineWidth',2);
set(gcf,'Color',[1 1 1]);
%end plotTransform_sub()

function [out]= makeTransform_sub (Gain, Bias, Linear)
%Create lookup table
if Bias < 0.00001
    Bias = 0.00001;
end;
if Gain < 0.00001
    Gain = 0.00001;
end;
in = 0:1/255:1;
out = 0:1/255:1;
if Linear
    midpoint = Bias;
    if Gain == 1.0
        slope = 1;
    else
        deg = Gain * 90; %0= no contrast = horizontal line, 1=no gray = vertical line
        rad = deg*pi/180; %degrees to radians
        slope = tan(rad); %http://en.wikipedia.org/wiki/Slope
    end;
    for i = 1:256,
        v = (((in(i) -midpoint)* slope)+ 0.5);
        
        if v > 1
            v = 1;
        elseif v < 0
                v = 0;
        end;
        out(i) = v;    
        %out(i) = in(i) * 2;
    end; %for i: each possible intensity
else %if not linear than non-linear
    %http://dept-info.labri.fr/~schlick/DOC/gem2.html
    %http://dept-info.labri.fr/~schlick/publi.html
    % "Fast Alternatives to Perlin's Bias and Gain Functions"
    %   Christophe Schlick: Graphics Gems IV, p379-382, April 1994 
    for i = 1:256,
        lT = in(i);
        %apply bias
        lT = (lT/((1/Bias-2)*(1-lT)+1)) ;
        %next apply gain
        gainX = 1.0-Gain;
        if lT <= 0.5 
            lG = (lT/((1/gainX-2)*(1-2*lT)+1));
        else
            lG = (( (1/gainX-2)*(1-2*lT)-lT ) / ( (1/gainX-2)*(1-2*lT)-1 ) );
        end;
        if lT == 0
            lG = 0;
        else
            lG = lG / lT;
        end;
        lV = (lT*lG);
        if lV > 1
            lV = 1;
        end;
        if lV < 0
            lV = 0;
        end;
        out(i) = lV;
        
    end; %for i: each possible intensity
end; %nonlinear 
%end makeTransform_sub()

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

function Im=RGB_YUV_sub(Y,U,V)
% This program transform YUV layers to RGB Layers in ome matrix 'Im'....
%  By   Mohammed Mustafa Siddeq
%  Date 25/7/2010
G=((Y-(U+V)/4));
R=U+G;
B=V+G;
Im(:,:,1)=R; Im(:,:,2)=G; Im(:,:,3)=B; 
%imshow(uint8(Im));
%end RGB_YUV_sub()

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
