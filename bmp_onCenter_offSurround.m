function bmp_onCenter_offSurround
%emulate response of on-center off-surround retinal ganglion cell
%http://agamtyagi.blogspot.com/2013/02/matlab-code-for-famous-mexican-hat.html  
kExtent = 6.5; %zeros at pi*2, 6.2832
kAmplitude = 20;
[x,y] = meshgrid(-kExtent:0.25:kExtent);
r = sqrt(x.^2+y.^2)+eps;
z = sin(r)./r;
z(r > (2*pi) & z > 0) = 0;
z = z * kAmplitude;
h = surf(z); %shading flat

set(h, 'edgecolor',[0 0 0], 'FaceLighting','phong');
set(h, 'FaceColor',[1 1 1], 'FaceAlpha',1, 'EdgeAlpha', 1);
set(gcf, 'color', [1 1 1])
axis off
axis equal