function bmp_sinewave_2d
%creates a 2D grayscale sinewave grating
%Clones Elliot Freeman code without requiring any toolboxes
%  http://www.icn.ucl.ac.uk/courses/MATLAB-Tutorials/Elliot_Freeman/html/gabor_tutorial.html

%parameters
contrast = 1.0; %contrast amplitude 0..1
imSize = 200;                           % image size: n X n
lamda = 20;                             % wavelength (number of pixels per cycle)
theta = 15;                              % grating orientation
sigma = 30;                             % gaussian standard deviation in pixels
phase = .25;                            % phase (0 -> 1)
trim = .005;                             % trim off gaussian values smaller than this
radial = false;                         %radial sine wave
%linear ramp
X = 1:imSize;                           % X is a vector from 1 to imageSize
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
freq = imSize/lamda;                    % compute frequency from wavelength
Xf = X0 * freq * 2*pi;                  % convert X to radians: 0 -> ( 2*pi * frequency)
sinX = sin(Xf) ;                        % make new sinewave
phaseRad = (phase * 2* pi);             % convert to radians: 0 -> 2*pi
[Xm Ym] = meshgrid(X0, X0);             % 2D matrices
Xf = Xm * freq * 2*pi;
grating = sin( Xf + phaseRad);          % make 2D sinewave
imagesc( grating, [-1 1] );             % display
colormap gray(256);                     % use gray colormap (0: black, 1: white)
set(0,'defaultaxesposition',[0 0 1 1])% display nicely without borders
%set(gcf, 'menu', 'none', 'Color',[.5 .5 .5]); % without background
%hFig = figure(1);
set(gcf, 'Position', [30 30 imSize imSize])
set(gca,'pos', [0 0 1 1]); 
axis off; axis image;                   % use gray colormap
if ~radial, return; end;%next bits for radial gradient
%make gaussian mask
s = sigma / imSize;                     % gaussian width as fraction of imageSize
%%Xg = exp( -( ( (X0.^2) ) ./ (2* s^2) ));% formula for 1D gaussian
%%Xg = normpdf(X0, 0, (20/imSize)); Xg = Xg/max(Xg);  % alternative using normalized probability function (stats toolbox)
%%plot(Xg)
gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
%imagesc( gauss, [-1 1] );                        % display
%axis off; axis image;     % use gray colormap
gauss = gauss * contrast; %contrast
gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)
gabor = grating .* gauss;                % use .* dot-product
imagesc( gabor, [-1 1] );                        % display
axis off; axis image;                    % use gray colormap
axis image; axis off; colormap gray(256);
%end bmp_sinewave_2d() main function - subfunctions follow

function pdf = normpdf (x, m, s)
%http://www.dynare.org/dynare-matlab-m2html/matlab/missing/stats/normpdf.html
if (nargin ~= 1 && nargin ~= 3)
 error('normpdf: you must give one or three arguments');
end
if (nargin == 1)
 m = 0;
 s = 1;
end
if (~isscalar (m) || ~isscalar (s))
 [retval, x, m, s] = common_size (x, m, s);
 if (retval > 0)
     error ('normpdf: x, m and s must be of common size or scalars');
 end
end
sz = size (x);
pdf = zeros (sz);
if (isscalar (m) && isscalar (s))
 if (find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf)))
     pdf = NaN * ones (sz);
 else
     pdf = stdnormal_pdf ((x - m) ./ s) ./ s;
 end
else
 k = find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf));
 if (any (k))
     pdf(k) = NaN;
 end
 k = find (~isinf (m) & ~isnan (m) & (s >= 0) & (s < Inf));
 if (any (k))
     pdf(k) = stdnormal_pdf ((x(k) - m(k)) ./ s(k)) ./ s(k);
 end
end
pdf((s == 0) & (x == m)) = Inf;
pdf((s == 0) & ((x < m) | (x > m))) = 0;
%end normpdf()

function pdf = stdnormal_pdf (x)
 %http://www.dynare.org/dynare-matlab-m2html/matlab/missing/stats/stdnormal_pdf.html
 if (nargin ~= 1)
     error('stdnormal_pdf: you should provide one argument');
 end
 
 sz = size(x);
 pdf = zeros (sz);
 
 k = find (isnan (x));
 if (any (k))
     pdf(k) = NaN;
 end
 
 k = find (~isinf (x));
 if (any (k))
     pdf (k) = (2 * pi)^(- 1/2) * exp (- x(k) .^ 2 / 2);
 end
%end stdnormal_pdf
