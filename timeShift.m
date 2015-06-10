function globalSignal = timeShift (restname, t1name, TRsec)
%Compute timeshift in each voxel
% restname: base name of resting state data
% t1name: base name of T1 scan: must be warped, segmented wc1 and wc2 images
% TRsec: repeat time for MRI data
%Example
% ts = timeShift('con001_a_rest.nii', 'con001_a_T1.nii', 2.0);
% plot(ts)
if ~exist('restname','var')
    restname = spm_select(1,'image','Select base 4D resting-state');
end
if ~exist('t1name','var')
    t1name = spm_select(1,'image','Select base T1-image');
end
if ~exist('TRsec','var') 
    TRsec = str2double(cell2mat(inputdlg('Repitition time (TR, sec)?:', 'Timing', 1,{'2'})));
end
gm = prefixSub ('wc1', t1name); %gray matter matter
rs = prefixSub ('fdwa', restname); %normalized resting state
globalSignal = timecourse(rs, gm); %global timecourse
wm = prefixSub ('wc2', t1name); %white matter map
mapShift(gm, wm, rs, globalSignal, TRsec);
%end timeShift()

function mapShift(gm, wm, rs, globalSignal, TRsec)
%1: make a masking image
kThresh = 0.25; %plot time shift for 
[hdr1, img1] = read_volsSub (gm); 
[hdr2, img2] = read_volsSub (wm);  %#ok<ASGLU>
imgTS = ((img1+ img2) > kThresh);
% make floating point
hdrTS = hdr1;
hdrTS.fname = prefixSub ('ts', rs);
hdrTS.dt(1) = 16; %make 32 bit real
hdrTS.private.dat.dtype = 'FLOAT32-LE';
hdrTS.private.dat.scl_slope = 1;
hdrTS.private.dat.scl_inter = 0;
hdrTS.pinfo = [1;0;352]; %slope=1, intercept = 0
%2: load resting state 
[hdr4d, img4d] = read_volsSub (rs); %#ok<ASGLU>
if numel(globalSignal) ~= size(img4d,4), error('Something is wrong'); end;
img4d = reshape(img4d, hdrTS.dim(1)*hdrTS.dim(2)*hdrTS.dim(3), size(img4d,4))'; %4D->2D
globalSignal = normColSub(globalSignal); %normalize signal to range 0..1
img4d = normColSub(img4d); %normalize each voxel to range 0..1
intercepts = regressSub(globalSignal, img4d); %compute intercepts
intercepts = reshape(intercepts, hdrTS.dim(1), hdrTS.dim(2), hdrTS.dim(3));%1D->3D
intercepts(imgTS == 0) = 0; %mask image
intercepts = intercepts * TRsec; %convert units to seconds
spm_write_vol(hdrTS,intercepts); %save data
%end mapShift()

function x = normColSub(x)
%normalize each column for range 0..1
% x = [1 4 3 0; 2 6 2 5; 3 10 2 2] -> x = [0 0 1 0; 0.5 0.333 0 1; 1 1 0 0.4]
x = bsxfun(@minus,x,min(x,[],1)); %translate so minimum = 0
x = bsxfun(@rdivide,x,max(x,[],1)); %scale so range is 1
%end normColSub()

function intercepts = regressSub(model, observed)
% gs = [1 2 3 4 5]';
% vox = [1.5 2.5 3.5 4.5 5.5; 1 2 3 4 5; 0 1 2 3 4]';
% regressSub(gs, vox)
X = [ model ones(length(model),1)]; %2 columns: slope and intercept
pXX = pinv(X)*pinv(X)'; % = pinv(X'*X), which is reusable, because
pX  = pXX * X';  % pinv(P*X) = pinv(X'*P'*P*X)*X'*P' = pXX * (P*X)'
b = pX * observed; %compute slopes and intercepts
intercepts = b(2,:);
%slopes = b(1,:); %<- unused, since we normalize values
%end regressSub

function signal = timecourse(img4d, imgMask)
%return mean signal in mask
[hdr4d, img4d] = read_volsSub (img4d); %#ok<ASGLU>
[hdrM, imgM] = read_volsSub (imgMask); %#ok<ASGLU>
if (size(imgM,1) ~= size(img4d,1)) || (size(imgM,2) ~= size(img4d,2)) || (size(imgM,3) ~= size(img4d,3))
    error('Image and mask must have same dimensions');
end
vols = size(img4d,4);
signal = zeros(vols, 1);
imgM(imgM == 0) = nan;
for v = 1 : vols
    i =img4d(:,:,:,v);
    i = i .* imgM;
    i = i(:);
    i = i(isfinite(i));
    signal(v) = mean(i);
end
%end timecourse

function [hdr, img] = read_volsSub (fnm)
[fnm, isGz] = unGzSub (fnm); %convert FSL .nii.gz to .nii
hdr = spm_vol(fnm); %load header data
img = spm_read_vols(hdr); %load image data
if (isGz), delete(fnm); end; %remove .nii if we have .nii.gz
%end read_volsSub()

function [fnm, isGz] = unGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm);
isGz = false;
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));  
    isGz = true;
elseif strcmpi(ext,'.voi') %.voi -> 
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
    isGz = true;
end;  
%end unGzSub()

function nam = prefixSub (pre, nam)
[p, n, x] = spm_fileparts(nam);
nam = fullfile(p, [pre, n, x]);
%end prefixSub()
