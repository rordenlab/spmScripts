function nii_cleanup5(gi,wi,ci, bi,si,ri);
%clean segmentation partitions for 5 tissue types, optionally create image for volume rendering
%  gi: gray matter image filename
%  wi: white matter image filename
%  ci: CSF matter image filename
%  bi: bone matter image filename [optional]
%  si: non-brain soft tissue image filename [optional]
%  ri: raw (T1-weighted) image from individual [optional]
%
%Example
%  nii_cleanup5('c1head_1.nii','c2head_1.nii','c3head_1.nii','c4head_1.nii','c5head_1.nii','head_1.nii')
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
if nargin <4 %no bone
 bi = spm_select(1,'image','Select bone volume');
end;
if nargin <5 %no non-brain soft tissue
 si = spm_select(1,'image','Select non-brain soft tissue volume');
end;

if ischar(gi), gi = spm_vol(gi); end;
if ischar(wi), wi = spm_vol(wi); end;
if ischar(ci), ci = spm_vol(ci); end;

nii_extract(gi.fname,0, wi.fname);
nii_extract(wi.fname);
nii_extract(ci.fname, 2500);%newborn eyeball 16.5mm 2352mm^3 [Riordan-Eva and Whitcher, 2008] average human eyeball 18mm diameter = 3000mm^3

[pth,nam,ext]=fileparts(gi.fname);
gi.fname = fullfile(pth,['e',  nam, ext]);
[pth,nam,ext]=fileparts(wi.fname);
wi.fname = fullfile(pth,['e',  nam, ext]);
[pth,nam,ext]=fileparts(ci.fname);
ci.fname = fullfile(pth,['e',  nam, ext]);

g = spm_read_vols(gi)*255;
w = spm_read_vols(wi)*255;
c = spm_read_vols(ci)*255;

if ischar(bi), bi = spm_vol(bi); end;
if ischar(si), si = spm_vol(si); end;
o = spm_read_vols(bi)*255;
s = spm_read_vols(si)*255;
o = o+s;


level = 2; %1=basic, 2=thorough
    
b    = w;

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
% for i=1:size(b,3),
%         gp       = double(g(:,:,i))/255;
%         wp       = double(w(:,:,i))/255;
%         cp       = double(c(:,:,i))/255;
%         op       = double(o(:,:,i))/255;
%         bp       = double(b(:,:,i))/255;
%         bp       = ((bp>th).*(wp+gp))>th;
%         g(:,:,i) = uint8(round(255*gp.*bp./(gp+wp+cp+eps)));
%         w(:,:,i) = uint8(round(255*wp.*bp./(gp+wp+cp+eps)));
%         c(:,:,i) = uint8(round(255*(cp.*bp./(gp+wp+cp+eps)+cp.*(1-bp))));
% end;
 for i=1:size(b,3),
         gp       = double(g(:,:,i))/255;
         wp       = double(w(:,:,i))/255;
         cp       = double(c(:,:,i))/255;
         op       = double(o(:,:,i))/255;
         bp       = double(b(:,:,i))/255;
         bp       = ((bp>th).*(wp+gp))>th;
         g(:,:,i) = uint8(round(255*gp.*bp./(gp+wp+cp+op+eps)));
         w(:,:,i) = uint8(round(255*wp.*bp./(gp+wp+cp+op+eps)));
         c(:,:,i) = uint8(round(255*(cp.*bp./(gp+wp+cp+op+eps)+cp.*(1-bp))));
 end;
spm_progress_bar('Clear');

spm_write_vol(gi,g/255);
spm_write_vol(wi,w/255);
spm_write_vol(ci,c/255);


if nargin <5 %no bone and soft tissue maps
 return;
end;

clear gp wp bp b 

if ischar(bi), bi = spm_vol(bi); end;
if ischar(si), si = spm_vol(si); end;

nii_extract(bi.fname);
nii_extract(si.fname);

[pth,nam,ext]=fileparts(bi.fname);
bi.fname = fullfile(pth,['e',  nam, ext]);
[pth,nam,ext]=fileparts(si.fname);
si.fname = fullfile(pth,['e',  nam, ext]);

c = g+w+c; % <-gray/white/CSF total...
g = spm_read_vols(bi)*255;
w = spm_read_vols(si)*255;

b    = w;

spm_progress_bar('Init',niter,'Extracting Bone/SoftTissue','Iterations completed');
for j=1:niter,
        if j>2, th=th1; else th=0.6; end; % Dilate after two its of erosion.
        for i=1:size(b,3),
                wp = double(w(:,:,i));
                bp = double(b(:,:,i))/255;
                bp = (bp>th).*(wp);
                b(:,:,i) = uint8(round(bp));
        end;
        spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
        spm_progress_bar('Set',j);
end;


th = 0.05;
for i=1:size(b,3),
        wp       = double(w(:,:,i))/255;
        bp       = double(b(:,:,i))/255;
        bp       = ((bp>th).*(wp))>th;
        w(:,:,i) = uint8(round(255*wp.*bp));
end;

spm_progress_bar('Clear');

spm_write_vol(bi,g/255);
spm_write_vol(si,w/255);

if nargin < 6 %no raw image
 return;
end;

%next lines make a brain render...
g = spm_read_vols(gi);
w = spm_read_vols(wi);
w = g+w;
mx=max(w(:));
w=w/mx;
if ischar(ri), ri = spm_vol(ri); end;
r = spm_read_vols(ri);
r=w.*r;
[pth,nam,ext]=fileparts(ri.fname);
ri.fname = fullfile(pth,['render',  nam, ext]);
spm_write_vol(ri,r);

%NEXT LINES CREATE BINARY BRAIN MASK
%r=round(w);
%ri.fname = fullfile(pth,['s',  nam, ext]);
%spm_write_vol(ri,r);
