function nii_cleanup3(gi,wi,ci, bi,si);
%clean segmentation partitions
%  gi: gray matter image filename
%  wi: white matter image filename
%  ci: CSF matter image filename
%
%Example
%  nii_cleanup('c1head_2.nii','c2head_2.nii','c3head_2.nii',2)
%
%Adapted John Ashburner's from spm_preproc_write
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

if nargin <1 %no gray
 gi = spm_select(1,'image','Select gray matter volume');
end;
if nargin <2 %no white matter
 wi = spm_select(1,'image','Select white matter volume');
end;
if nargin <3 %no CSF
 ci = spm_select(1,'image','Select CSF volume');
end;

if ischar(gi), gi = spm_vol(gi); end;
if ischar(wi), wi = spm_vol(wi); end;
if ischar(ci), ci = spm_vol(ci); end;

g = spm_read_vols(gi)*255;
w = spm_read_vols(wi)*255;
c = spm_read_vols(ci)*255;

level = 2; %1=basic, 2=thorough
   
b    = w;
b(1) = w(1);

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level==2, th1 = 0.2; end;
% Erosions and conditional dilations
%--------------------------------------------------------------------------
niter = 32;
   
spm_progress_bar('Init',niter,'Extracting Brain','Iterations completed');
for j=1:niter,
        if j>2, th=th1; else th=0.6; end; % Dilate after two its of erosion.
        for i=1:size(b,3),
                gp = double(g(:,:,i));
                wp = double(w(:,:,i));
                bp = double(b(:,:,i))/255;
                bp = (bp>th).*(wp+gp);
                b(:,:,i) = uint8(round(bp));
        end;
        spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
        spm_progress_bar('Set',j);
end;

th = 0.05;
for i=1:size(b,3),
        gp       = double(g(:,:,i))/255;
        wp       = double(w(:,:,i))/255;
        cp       = double(c(:,:,i))/255;
        bp       = double(b(:,:,i))/255;
        bp       = ((bp>th).*(wp+gp))>th;
        g(:,:,i) = uint8(round(255*gp.*bp./(gp+wp+cp+eps)));
        w(:,:,i) = uint8(round(255*wp.*bp./(gp+wp+cp+eps)));
        c(:,:,i) = uint8(round(255*(cp.*bp./(gp+wp+cp+eps)+cp.*(1-bp))));
end;
spm_progress_bar('Clear');

[pth,nam,ext]=fileparts(gi.fname);
gi.fname = fullfile(pth,['e',  nam, ext]);
spm_write_vol(gi,g/255);

[pth,nam,ext]=fileparts(wi.fname);
wi.fname = fullfile(pth,['e',  nam, ext]);
spm_write_vol(wi,w/255);

[pth,nam,ext]=fileparts(ci.fname);
ci.fname = fullfile(pth,['e',  nam, ext]);
spm_write_vol(ci,c/255);
