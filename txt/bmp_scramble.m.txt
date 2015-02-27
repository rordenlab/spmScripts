function bmp_scramble (Filename)
%Creates phase spectrum scrambling bitmap with 's' prefix
%  Filename: name of bitmap [optional]
%Example
% bmp_scramble('dog.png');
%Adds a user interface wrapper for Nicolaas Prins' code
% http://visionscience.com/pipermail/visionlist/2007/002181.html
%modified by Chris Rorden so image processing toolbox is not required

if (nargin < 1)  
   [files,pth] = uigetfile({'*.bmp;*.jpg;*.png;*.tiff;';'*.*'},'Select the Image[s]', 'MultiSelect', 'on'); 
   files = cellstr(files); %make cellstr regardless of whether user selects single or multiple images
else
    [pth,nam, ext] = fileparts(Filename);
    files = cellstr([nam, ext]); 
end;

%if (license('checkout', 'signal_toolbox')) == 0 
%    error('bmp_scramble error: You will need the image processing toolbox!')
%end;

for i=1:size(files,2)
    nam = strvcat(deblank(files(:,i))); %#ok<REMFF1>
    Inname = fullfile(pth, nam);
    Outname = fullfile(pth, ['s' nam ]);
    %Im =  mat2gray(double(imread(Inname)));% <- requires signal processing
    Im = Imread_mat2gray_sub(Inname); 
    tic;
    %rand('twister',76599) % <- set random seed to replicate results
    ImSize = size(Im);
    RandomPhase = angle(fft2(rand(ImSize(1), ImSize(2))));
    %determine layers: grayscale =1, red/green/blue =3
    if length(ImSize) == 2
        n = 1; %only one layer - e.g. grayscale image
    else
        n = ImSize(3); %multiple layers, e.g. color image
    end;
    %generate random phase structure
    ImFourier = Im; %preallocate Fourier
    Amp = Im; %preallocate Amplitude
    Phase = Im; %preallocate Phase
    ImScrambled = Im;
    for layer = 1:n
        ImFourier(:,:,layer) = fft2(Im(:,:,layer));       
        %Fast-Fourier transform
        Amp(:,:,layer) = abs(ImFourier(:,:,layer));       
        %amplitude spectrum
        Phase(:,:,layer) = angle(ImFourier(:,:,layer));   
        %phase spectrum
        Phase(:,:,layer) = Phase(:,:,layer) + RandomPhase;
        %add random phase to original phase
        ImScrambled(:,:,layer) = ifft2(Amp(:,:,layer).*exp(sqrt(-1)*(Phase(:,:,layer))));   
        %combine Amp and Phase then perform inverse Fourier
    end
    ImScrambled = real(ImScrambled); %get rid of imaginary part in image (due to rounding error)
    toc
    imwrite(ImScrambled,Outname);
    %imshow(ImScrambled) %display result: requires image processing toolbox
    clear Amp Phase Im ImFourier ImScrambled %free these incase subsequent images have different sizes
end;

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
    result = Im; %preallocate image
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
%end Imread_mat2gray_sub

