function t1Bet = extractSub(thresh, t1, probMap)   
%subroutine to extract brain from surrounding scalp
% t1: image extracted, requires corresponding 'c1','c2' tissue maps
fprintf('Brain extraction of %s\n', t1);
[pth,nam,ext] = spm_fileparts(t1);
c1 = fullfile(pth,['c1',nam,ext]);
c2 = fullfile(pth,['c2',nam,ext]);
if ~exist(c1,'file') || ~exist(c2,'file')
    error('Unable to find tissue maps %s %s', c1,c2);
end
%load headers
mi = spm_vol(t1);%bias corrected T1
gi = spm_vol(c1);%Gray Matter map
wi = spm_vol(c2);%White Matter map
%load images
m = spm_read_vols(mi);
g = spm_read_vols(gi);
w = spm_read_vols(wi);
%w = g+w;
max(w(:))
if thresh <= 0
    if ~exist('probMap','var') || ~probMap
        m=m.*w;
    end;
else
    mask= zeros(size(m));
    mask(w >= thresh) = 255;
    maskIn = mask;
    spm_smooth(maskIn,mask,1); %feather the edges
    mask = mask / 255;
    if ~exist('probMap','var') || ~probMap
        m=m.*mask;
    else
        m = mask;
    end;
    
end;
if exist('probMap','var') && probMap
    m = m * 1000;
    mi.fname = fullfile(pth,['p',  nam, ext]);
else
    mi.fname = fullfile(pth,['b',  nam, ext]);

end
mi.dt(1) = 4; %16-bit precision more than sufficient uint8=2; int16=4; int32=8; float32=16; float64=64
spm_write_vol(mi,m);
t1Bet = mi.fname;
%end extractSub()
