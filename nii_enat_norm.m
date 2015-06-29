function nii_enat_norm(T1,lesion,T2, UseXTemplate, vox, bb, DeleteIntermediateImages, ssthresh, autoOrigin)
%Perform enantiomorphic normalization using SPM12
% see Nachev et al. (2008) http://www.ncbi.nlm.nih.gov/pubmed/18023365
%  T1: filename of T1 image
%  Lesion: filename of lesion map 
%  T2: (optional) filename of image used to draw lesion, if '' then lesion drawn on T1
%  UseXTemplate: if false (default) standard SPM template is used, else special template 
%Examples
% nii_enat_norm('T1_LM1054.nii','LS_LM1054.nii','') %lesion drawn on T1 scan
% nii_enat_norm('T1_LM1054.nii','LS_LM1054.nii',''T2_LM1054.nii') %lesion drawn on T2 scan
% nii_norm12('MB_T1.nii','',''); %no lesion - control participant
% nii_norm12 %use graphical interface

%STEP 0: check inputs
isSPM12orNewerSub;spm fmri
T1param = exist('T1','var'); %did the user provide a T1 scan
if ~T1param, T1 = spm_select(1,'image','Select T1 images'); end;
if isempty(T1), return; end;
if ~exist('lesion','var') && ~T1param, lesion = spm_select(1,'image','Optional: select lesion map'); end;
if ~exist('lesion','var'), lesion = ''; end;
if ~isempty(lesion) && ~exist('T2','var'), T2 = spm_select(1,'image','Optional: Select image used to draw lesion (if not T1)'); end;
if ~exist('T2','var'), T2 = ''; end;
if ~exist('UseXTemplate','var'),UseXTemplate = 0; end;
if ~exist('vox','var'), vox = [1 1 1]; end;
if ~exist('bb','var'), bb = [-78 -112 -70; 78 76 85]; end;
if ~exist('DeleteIntermediateImages','var'), DeleteIntermediateImages = true; end;
if ~exist('ssthresh','var'), ssthresh = 0.005; end; %with SPM12, better GM, so threshold of 1%
if ~exist('autoOrigin','var')
   %ButtonName = questdlg('Automatic origin detection?','Preferences', 'Yes', 'No', 'No');
   %autoOrigin = strcmpi(ButtonName,'Yes');
   autoOrigin = false;
end
T1 = stripVolSub(T1); lesion = stripVolSub(lesion); T2 = stripVolSub(T2);
if isDoneSub(T1), fprintf('Already done: skipping normalization of %s\n',T1); return; end;
[T1,lesion,T2] = checkDimsSub(T1, lesion, T2); %check alignment
%0: rough estimate for origin and alignment
if autoOrigin
    setOriginSub({T1, T2, lesion}, 1); 
end
%1: align lesion/t2 to match T1
[rT2, rlesion] = coregEstWriteSub(T1,T2,lesion); %#ok<ASGLU>

rlesion = smoothSub(rlesion, 3); 
%2: make image without lesion
eT1 = entiamorphicSub (T1, rlesion);
%[eT1, erT2] = entiamorphicSub (T1, rlesion, rT2); %for multichannel
%3: new-segment image
newSegSub(eT1,'', UseXTemplate);
%newSegSub(eT1, erT2, UseXTemplate); %for multichannel
%4: create 'b' (brain extracted) image without scalp signal
bT1 = extractSub(ssthresh, T1, prefixSub('c2', eT1), prefixSub('c1', eT1));
%5: warp render image to standard space
rT1 = newSegWriteSub(eT1, bT1, [0.9 0.9 0.9]); %#ok<NASGU>
%6: warp lesion to standard space
wrlesion = newSegWriteSub(eT1, rlesion, vox, bb, true); %#ok<NASGU>
wT1 = newSegWriteSub(eT1, T1, vox, bb); %#ok<NASGU>
if DeleteIntermediateImages, deleteSub(T1); end;
%end nii_enat_norm()

%--- local functions follow 
%function img = smoothSub(img, FWHM)
%[pth,nam,ext] = spm_fileparts(img);
%smth = fullfile(pth, ['s' nam ext]);
%spm_smooth(img, smth, FWHM, 0);  
%img = smth;
%end smoothSub()

function isDone = isDoneSub(T1)
isDone = false;
[pth,nam,ext] = fileparts(T1);
b = fullfile(pth,['b', nam, ext]); %brain extracted image
if ~exist(b,'file'), return; end; 
defname = fullfile(pth,['y_' nam ext]); %control normalization
edefname = fullfile(pth,['y_e' nam ext]); %patient normalization
if exist(defname,'file') || exist(edefname,'file'), isDone = true; end;
%end isDoneSub()

function img = smoothSub(img, FWHM)
if isempty(img), return; end;
hdr = spm_vol(img);
im = spm_read_vols(hdr);
if (spm_type(hdr.dt,'intt')) %integer data
	mn = min(im(:));
	range = max(im(:)) - mn;
	if range < 10 && range > 0
		im = (im - mn) * 255/range;
	end
	hdr.pinfo(1) = 1; %slope
	hdr.pinfo(2) = 0; %intercept
end
smoothFWHMmm = [FWHM FWHM FWHM];
VOX = sqrt(sum(hdr.mat(1:3,1:3).^2));
smoothFWHMvox = smoothFWHMmm/VOX; %for 3D arrays the FWHM is specified in voxels 
presmooth = im+0; %+0 forces new matrix
spm_smooth(presmooth,im,smoothFWHMvox,0);
[pth,nam,ext] = spm_fileparts(img);
img = fullfile(pth, ['s' nam ext]);
hdr.fname = img;
spm_write_vol (hdr, squeeze (im ));

function img = stripVolSub(img)
%strip volume from lesion name, 'img.nii,1' -> 'img.nii'
if isempty(img), return; end;
[n,m,x] = spm_fileparts(img); %we ignore the volume
img = fullfile(n, [m, x]);
%end stripVolSub()

function [T1,lesion,T2] = checkDimsSub(T1, lesion, T2)
if ~exist(T1,'file'), error('T1 image required %s', T1); end;
if ~exist('lesion','var'), return; end;
if ~exist(lesion,'file'), error('Lesion image not found %s', lesion); end;
hdrT1 = spm_vol(T1);
hdrLS = spm_vol(lesion);
mmLS = (hdrLS.mat * [0 0 0 1]'); %vox2mm,  [0 0 0 1; 0 0 1 1]'
mmT1 = (hdrT1.mat * [0 0 0 1]');
dxT1 = sqrt(sum((mmT1(1:3)-mmLS(1:3)).^2)); %error between T1 and lesion
if exist('T2','var') && ~isempty(T2)
   if ~exist(T2,'file'), error('T2 image not found %s', T2); end;
   hdrT2 = spm_vol(T2); 
   %we compute distance differently for T2/Lesion as these will be resliced...
   [mnT2, mxT2] = bbSub(hdrT2); %range of T2 bounding box
   [mnLS, mxLS] = bbSub(hdrLS); %range of Lesion bounding box
   dxMn = sqrt(sum((mnT2(1:3)-mnLS(1:3)).^2)); %error between T2 and lesion
   dxMx = sqrt(sum((mxT2(1:3)-mxLS(1:3)).^2)); %error between T2 and lesion
   if (dxMn > 1) || (dxMx > 1)
        if (dxT1 < 0.25)
            T2 = '';
            fprintf('WARNING: T2 dimensions do not match lesion - ASSUME lesion drawn on T1');
        else
            fprintf('WARNING: Neither T2 nor T1 aligned to lesion.');
        end 
   end
   return;
end %if T2 is present
if ~all(hdrT1.dim == hdrLS.dim)
    error('WARNING: T1 dimensions do not match lesion %s %s',T1, lesion);    
end
if (dxT1 > 0.25)
    fprintf('WARNING: T1 poorly aligned to lesion.');
end
%end checkDimsSub()

function [mn, mx] = bbSub(hdr) %return range for image bounding box
d = hdr.dim(1:3);
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
tc = hdr.mat(1:3,1:4)*c;
% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
%end bbSub()

function deleteSub(T1)
deletePrefixSub ('LR', T1)
deletePrefixSub ('rLR', T1)
%deleteSub

function deletePrefixSub (pre, nam)
nam = prefixSub (pre, nam);
if exist(nam, 'file'), delete(nam); end;
%end deletePrefixSub()

function isSPM12orNewerSub
%check that SPM is installed and is at least release 6225
if exist('spm','file') ~= 2, error('Please install SPM12 or later'); end;
[v,r] = spm('Ver','',1); r = str2double(r); %#ok<ASGLU>
if r < 6225, error('Please update your copy of SPM'); end;
%end isSPM12orNewer()

function t1Bet = extractSub(thresh, t1, c1, c2, c3)   
%subroutine to extract brain from surrounding scalp
% t1: anatomical scan to be extracted
% c1: gray matter map
% c2: white matter map
% c3: [optional] spinal fluid map
fprintf('Brain extraction of %s\n', t1);
[pth,nam,ext] = spm_fileparts(t1);
%load headers
mi = spm_vol(t1);%bias corrected T1
gi = spm_vol(c1);%Gray Matter map
wi = spm_vol(c2);%White Matter map
%load images
m = spm_read_vols(mi);
g = spm_read_vols(gi);
w = spm_read_vols(wi);
if nargin > 4 && ~isempty(c3)
   ci = spm_vol(c3);%CSF map
   c = spm_read_vols(ci);
   w = c+w; 
end;
w = g+w;
if thresh <= 0
    m=m.*w;
else
    mask= zeros(size(m));
    for px=1:length(w(:)),
      if w(px) >= thresh
        mask(px) = 255;
      end;
    end;
    spm_smooth(mask,mask,1); %feather the edges
    mask = mask / 255;
    m=m.*mask;
end;
mi.fname = fullfile(pth,['b',  nam, ext]);
mi.dt(1) = 4; %16-bit precision more than sufficient uint8=2; int16=4; int32=8; float32=16; float64=64
spm_write_vol(mi,m);
t1Bet = mi.fname;
%end extractSub()

function  targetname = newSegWriteSub(t1name, targetname, vox, bb, binarize)
%reslice img using pre-existing new-segmentation deformation field
if isempty(targetname) || isempty(t1name), return; end;
if ~exist('bb','var'), bb = [-78 -112 -50; 78 76 85]; end;
if ~exist('vox','var'), vox =[2 2 2]; end;
[pth,nam,ext, vol] = spm_fileparts(t1name); %#ok<NASGU>
defname = fullfile(pth,['y_' nam ext]);
if ~exist(defname,'file')
    error('Unable to find new-segment deformation image %s',defname);
end
fprintf('Warping %s based on NewSegment of %s\n', targetname, t1name);
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {defname};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {targetname};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = bb;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vox;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1; %4; //trilinear avoids ringing
spm_jobman('run',matlabbatch);
targetname = prefixSub('w', targetname); 
if ~exist('binarize','var') || ~binarize, return; end;
hdr = spm_vol(targetname);
img = spm_read_vols(hdr);
mn = min(img(:));
mx = max(img(:));
thresh = ((mx-mn)*0.5) + mn;
spm_write_vol(hdr,+(img > thresh));
%end newSegWriteSub()

function newSegSub(t1, t2, UseXTemplate)
%apply new segment - return name of warping matrix
template = fullfile(spm('Dir'),'tpm','TPM.nii');
if nargin > 2 && UseXTemplate
    xtemplate = fullfile(spm('Dir'),'toolbox','Clinical','TPM4mm.nii');
    if exist(xtemplate,'file')
        template = xtemplate;
    else
        fprintf('WARNING: unable to find template named %s\n', xtemplate);
    end
end
if ~exist(template,'file')
    error('Unable to find template named %s',template);
end
fprintf('NewSegment of %s\n', t1);
matlabbatch{1}.spm.spatial.preproc.channel(1).vols = {t1};
matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];
if nargin > 1 && ~isempty(t2)
    matlabbatch{1}.spm.spatial.preproc.channel(2).vols = {t2};
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel(2).write = [0 0];
end;
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[template ',1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[template ',2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[template ',3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[template ',4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[template ',5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[template ',6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
spm_jobman('run',matlabbatch);
%end newSegSub()

function [intactT1, intactT2] = entiamorphicSub (T1, lesionNam, T2)
%Generates image suitable for Enantiomorphic normalization, see www.pubmed.com/18023365
% anatNam   : filename of anatomical scan
% lesionNam : filename of lesion map in register with anatomical
%returns name of new image with two 'intact' hemispheres
if ~exist('T1','var') %no files specified
    T1 = spm_select(1,'image','Select anatomical image');
end
if ~exist('lesionNam','var') %no files specified
    lesionNam = spm_select(1,'image','Select lesion image');
end
if isempty(lesionNam) 
    intactT1 = T1;
    return;
end
if (exist(T1,'file') == 0) || (exist(lesionNam,'file') == 0)
    error('%s unable to find files %s or %s',mfilename, T1, lesionNam);
end
if nargin < 3, T2 = ''; end;
fprintf('Using lesion %s to substitute %s\n', lesionNam, T1);
%create flipped image
T1lr = flipSub(T1);
T2lr = flipSub(T2);
[T1lr, T2lr] = coregEstWriteSub(T1, T1lr, T2lr); %reslice mirror
intactT2 = insertSub(T2, T2lr, lesionNam);
intactT1 = insertSub(T1, T1lr, lesionNam);
%end entiamorphicSub()

function namFilled = insertSub(nam, namLR, lesion)
%namLR donates voxels masked by lesion to image nam
if isempty(nam), namFilled =''; return; end;
hdrLesion = spm_vol(lesion); 
imgLesion = spm_read_vols(hdrLesion);
rdata = +(imgLesion > (max(imgLesion(:))/2)); %binarize raw lesion data, + converts logical to double
spm_smooth(rdata,imgLesion,4); %blur data
rdata = +(imgLesion > 0.05); %dilate: more than 5%
spm_smooth(rdata,imgLesion,8); %blur data
%now use lesion map to blend flipped and original image
hdr = spm_vol(nam); 
img = spm_read_vols(hdr);
hdr_flip = spm_vol(namLR); 
imgFlip = spm_read_vols(hdr_flip);
size(img)
size(imgLesion)
rdata = (img(:) .* (1.0-imgLesion(:)))+ (imgFlip(:) .* imgLesion(:));
rdata = reshape(rdata, size(img));
[pth, nam, ext] = spm_fileparts(hdr.fname);
hdr_flip.fname = fullfile(pth,['e' nam ext]);%image with lesion filled with intact hemisphere
spm_write_vol(hdr_flip,rdata); 
namFilled = hdr_flip.fname;
%insertSub()

function namLR = flipSub (nam)
if isempty(nam), namLR = ''; return; end;
hdr = spm_vol(nam); 
img = spm_read_vols(hdr);
[pth, nam, ext] = spm_fileparts(hdr.fname);
namLR = fullfile(pth, ['LR', nam, ext]);
hdr_flip = hdr;
hdr_flip.fname = namLR;
hdr_flip.mat = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * hdr_flip.mat;
spm_write_vol(hdr_flip,img); 
%end flipSub()

function [T2, lesion] = coregEstWriteSub(T1, T2, lesion)
%coregister T2 to match T1 image, apply to lesion
if isempty(T1) || isempty(T2), return; end;
fprintf('Coregistering %s to match %s\n',T2,T1);
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {T1};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {T2};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {lesion};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
T2 = prefixSub('r',T2);
if ~isempty(lesion), lesion = prefixSub('r',lesion); end;
%end coregEstSub()

function nam = prefixSub (pre, nam)
[p, n, x] = spm_fileparts(nam);
nam = fullfile(p, [pre, n, x]);
%end prefixSub()

function coivox = setOriginSub(vols, modality)
%Align images so that origin and alignment roughly match MNI space
%  vols : cell string of image name(s) - first image used for estimate, others yoked
%  modality : modality of first image 1=T1, 2=T2, 3=EPI
%Example
%  setOrigin('T1.nii',1); %align T1 scan
%  setOrigin({'T1s005.nii', 'fmriblocks009.nii'},1); %use T1 to align T1 and fMRI data
%  setOrigin %use graphical interface
%Chris Rorden 12/2014 (now supports SPM12)
if ~exist('vols','var') %no files specified
 vols = spm_select(inf,'image','Reset origin for selected image(s) (estimated from 1st)');
end
if ischar(vols)
    vols = cellstr(vols);
end
if ~exist('modality','var') %no files specified
 modality = 1;
 fprintf('%s Modality not specified, assuming T1\n', mfilename);
end
coivox = ones(4,1);
%extract filename 
[pth,nam,ext, ~] = spm_fileparts(deblank(vols{1}));
fname = fullfile(pth,[nam ext]); %strip volume label
%report if filename does not exist...
if (exist(fname, 'file') ~= 2) 
 	fprintf('%s error: unable to find image %s.\n',mfilename,fname);
	return;  
end;
hdr = spm_vol([fname,',1']); %load header 
img = spm_read_vols(hdr); %load image data
img = img - min(img(:));
img(isnan(img)) = 0;
%find center of mass in each dimension (total mass divided by weighted location of mass
% img = [1 2 1; 3 4 3];
sumTotal = sum(img(:));
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZ_mm = hdr.mat * coivox; %convert from voxels to millimeters
fprintf('%s center of brightness differs from current origin by %.0fx%.0fx%.0fmm in X Y Z dimensions\n',fname,XYZ_mm(1),XYZ_mm(2),XYZ_mm(3)); 
for v = 1:   numel(vols) 
    fname = deblank(vols{v});
    if ~isempty(fname)
        [pth,nam,ext, ~] = spm_fileparts(fname);
        fname = fullfile(pth,[nam ext]); 
        hdr = spm_vol([fname ',1']); %load header of first volume 
        fname = fullfile(pth,[nam '.mat']);
        if exist(fname,'file')
            destname = fullfile(pth,[nam '_old.mat']);
            copyfile(fname,destname);
            fprintf('%s is renaming %s to %s\n',mfilename,fname,destname);
        end
        hdr.mat(1,4) =  hdr.mat(1,4) - XYZ_mm(1);
        hdr.mat(2,4) =  hdr.mat(2,4) - XYZ_mm(2);
        hdr.mat(3,4) =  hdr.mat(3,4) - XYZ_mm(3);
        spm_create_vol(hdr);
        if exist(fname,'file')
            delete(fname);
        end
    end
end%for each volume
coregEstTemplateSub(vols, modality);
for v = 1:   numel(vols) 
    [pth, nam, ~, ~] = spm_fileparts(deblank(vols{v}));
    fname = fullfile(pth,[nam '.mat']);
    if exist(fname,'file')
        delete(fname);
    end
end %for each volume
%end setOriginSub()


function coregEstTemplateSub(vols, modality)
if modality == 2
   template = fullfile(spm('Dir'),'canonical','avg152T2.nii');
elseif modality == 3
    template  = fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii');
else
    template = fullfile(spm('Dir'),'canonical','avg152T1.nii');
end
if ~exist(template,'file')
    error('Unable to find template named %s\n', template);
end
if ischar(vols)
    vols = cellstr(vols);
end
vols(strcmp('',vols)) = []; %remove empty strings
%matlabbatch{1}.spm.spatial.coreg.estimate.ref = {'/Users/Shared/spm12/canonical/avg152PD.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {template};
%matlabbatch{1}.spm.spatial.coreg.estimate.source = {'/Users/rorden/Desktop/pre/bvisiblehuman.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[deblank(vols{1}),',1']};%{'/Users/rorden/Desktop/3D.nii,1'};
if  numel(vols) > 1
    matlabbatch{1}.spm.spatial.coreg.estimate.other = vols(2:end);% {''};
else
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
end
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
%end coregEstTemplateSub()