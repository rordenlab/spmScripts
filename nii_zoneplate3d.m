function nii_zoneplate3d(N)
%make Fresnel zone plates - useful for checking aliasing of resampling methods
% https://en.wikipedia.org/wiki/Zone_plate
% https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/35961/versions/2/previews/imzoneplate.m/index.html?access_key=
% https://www.mathworks.com/matlabcentral/fileexchange/35961-zone-plate-test-image

if ~exist('N','var')    
    N = 129;
end
hdr.fname = sprintf('zoneplate3d_%d.nii', N);  
%create header
hdr.dim = [N, N, N];
hdr.dt = [2 0];
hdr.pinfo = [0;0;0];
hdr.mat = [1 0 0 -N/2; 0 1 0 -N/2; 0 0 1 -N/2; 0 0 0 1];
hdr.descrip = 'https://en.wikipedia.org/wiki/Zone_plate';
%create image
img = zeros(N,N,N);
if rem(N,2) == 1
    x2 = (N-1)/2;
    x1 = -x2;
else
    x2 = N/2;
    x1 = -x2 + 1;
end
x = [x1:x2];
for z = 1:N
    img(:,:,z) = imzoneplate(N, x(z));  
end
%save image
spm_write_vol(hdr,img);



function I = imzoneplate(N, dx)
%imzoneplate Zone plate test pattern
%
%   SYNTAX
%
%   I = imzoneplate
%   I = imzoneplate(N)
%
%   DESCRIPTION
%
%   I = imzoneplate creates a 501-by-501 zone plate test image. This is a
%   radially symmetric pattern with low frequencies in the middle and high
%   frequencies near the edge.
%
%   I = imzoneplate(N) creates an N-by-N zone plate test image.
%
%   EXAMPLES
%
%   Create a test image with the default size (501-by-501).
%
%       I = imzoneplate;
%       imshow(I)
%
%   Create a smaller test image and plot its cross-section.
%
%       I = imzoneplate(151);
%       plot(I(76,:))
%
%   REFERENCE
%
%   Practical Handbook on Image Processing for Scientific Applications, by
%   Bernd Jhne, CRC Press, 1997. See equation 10.63:
%
%   g({\bf x}) = g_0 \sin\left(\frac{k_m|{\bf x}|^2}{2r_m}\right) 
%   \left[\frac{1}{2} \tanh\left(\frac{r_m-|{\bf x}|}{w}\right) + 
%   \frac{1}{2}\right]
%
%   In this equation, g takes on values in the range [-1,1]. imzoneplate
%   returns I = (g+1)/2, which takes on values in the range [0,1].
%
%   See also http://blogs.mathworks.com/steve/2011/07/19/jahne-test-pattern-take-3/

%   Copyright 2012 The MathWorks, Inc.
%   Steven L. Eddins

if nargin < 1
    N = 501;
end

if rem(N,2) == 1
    x2 = (N-1)/2;
    x1 = -x2;
else
    x2 = N/2;
    x1 = -x2 + 1;
end
[x,y] = meshgrid(x1:x2);
if ~exist('dx','var') || (dx == 0)
    r = hypot(x,y);
else
    r = zeros(N,N);
    dx2 = dx ^2;
    for i = 1 : numel(x)
        r(i) = sqrt(x(i)^2 + y(i)^2 + dx2);
    end 
end
%I = r; return
km = 0.7*pi;
rm = x2;
w = rm/10;
term1 = sin( (km * r.^2) / (2 * rm) );
term2 = 0.5*tanh((rm - r)/w) + 0.5;
g = term1 .* term2;

I = (g + 1)/2;

