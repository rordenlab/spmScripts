function nii_cleanup1(wi, level);
%clean segmentation partition
%  wi: image to cleanup
%Example:
%  nii_cleanup1('c1head_2.nii')
%Adapted from John Ashburner's spm_preproc_write.m

if nargin <1 %no gray
 wi = spm_select(1,'image','Select volume to cleanup');
end;
if nargin <2 %no gray
    level = 1; %1=basic, 2=thorough
end;

if ischar(wi), wi = spm_vol(wi); end;

w = spm_read_vols(wi)*255;
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

[pth,nam,ext]=fileparts(wi.fname);
wi.fname = fullfile(pth,['z',  nam, ext]);
spm_write_vol(wi,w/255);

