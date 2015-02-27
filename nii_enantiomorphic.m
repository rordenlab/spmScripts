function nii_enantiomorphic (anat,lesion)
%Enantiomorphic normalization, see Nachev et al (2008) http://www.ncbi.nlm.nih.gov/pubmed/18023365
% anat : filename(s) for anatomical scans
% lesion : filename(s) for lesion maps
%Chris Rorden 2014
% http://opensource.org/licenses/BSD-2-Clause
%Examples
% nii_enantiomorphic; %use GUI
% nii_enantiomorphic('T1.nii','lesion.nii');

fprintf('This is a minimal script for enantiomophic normalization: the clinical toolbox provides many enhancements.\n');
if exist('spm','file') ~= 2
    error('Please install SPM');
end
if ~exist('anat','var') %no files specified
    anat = spm_select(inf,'image','Select anatomical scan(s)');
end
if ~exist('lesion','var') %no files specified
    lesion = spm_select(inf,'image','Select lesion image(s) [same order]');
end
if (size(anat,1) < 1) || (size(anat,1) ~= size(lesion,1))
    error('Please specify images');
end
for i = 1: size(anat,1)
    fnm = deblank(anat(i,:));
    les = deblank(lesion(i,:));
    %make image with homologous healthy tissue inserted in place of lesion
    efnm = entiamorphicSub(fnm,les);
    %normalize
    nbatch{1}.spm.spatial.preproc.data = {efnm};
    nbatch{1}.spm.spatial.preproc.output.GM = [0 0 0];
    nbatch{1}.spm.spatial.preproc.output.WM = [0 0 0];
    nbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];
    nbatch{1}.spm.spatial.preproc.output.biascor = 0;
    nbatch{1}.spm.spatial.preproc.output.cleanup = 0;
    nbatch{1}.spm.spatial.preproc.opts.tpm = { fullfile(spm('Dir'),'tpm','grey.nii');fullfile(spm('Dir'),'tpm','white.nii');fullfile(spm('Dir'),'tpm','csf.nii') };
    nbatch{1}.spm.spatial.preproc.opts.ngaus = [2; 2; 2; 4];
    nbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
    nbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
    nbatch{1}.spm.spatial.preproc.opts.warpco = 25;
    nbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
    nbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
    nbatch{1}.spm.spatial.preproc.opts.samp = 3;
    nbatch{1}.spm.spatial.preproc.opts.msk = {''};
    spm_jobman('run',nbatch);
    %reslice anatomical and lesion images
    [pth,nam] = spm_fileparts(efnm);
    rbatch{1}.spm.spatial.normalise.write.subj.matname = {fullfile(pth,[ nam '_seg_sn.mat'])};
    rbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
	rbatch{1}.spm.spatial.normalise.write.roptions.interp = 1;
	rbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
	rbatch{1}.spm.spatial.normalise.write.subj.resample =  {fnm ,',1; ', les,',1'};
    rbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
    rbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50; 78 76 85];
    rbatch{1}.spm.spatial.normalise.write.roptions.vox = [2 2 2];
    spm_jobman('run',rbatch);
end

function intactImg = entiamorphicSub (anatImg, lesionImg)
%Generates image suitable for Enantiomorphic normalization, see www.pubmed.com/18023365
% anatImg   : filename of anatomical scan
% lesionImg : filename of lesion map in register with anatomical
%returns name of new image with two 'intact' hemispheres
if (exist(anatImg,'file') == 0) || (exist(lesionImg,'file') == 0)
    error('%s unable to find files %s or %s',mfilename, anatImg, lesionImg);
end
%create flipped image 
hdr = spm_vol(anatImg); 
img = spm_read_vols(hdr);
[pth, nam, ext] = spm_fileparts(hdr.fname);
fname_flip = fullfile(pth, ['LR', nam, ext]);
hdr_flip = hdr;
hdr_flip.fname = fname_flip;
hdr_flip.mat = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * hdr_flip.mat;
spm_write_vol(hdr_flip,img); 
%coregister data
hdr_flip = spm_vol(fname_flip); 
x  = spm_coreg(hdr_flip,hdr); 
%apply half of transform to find midline
x  = (x/2); 
M = spm_matrix(x);
MM = spm_get_space(fname_flip);
spm_get_space(fname_flip, M*MM); %reorient flip
M  = inv(spm_matrix(x)); 
MM = spm_get_space(hdr.fname);
spm_get_space(hdr.fname, M*MM); %#ok<MINV> %reorient original so midline is X=0
%reorient the lesion as well
MM = spm_get_space(lesionImg);
spm_get_space(lesionImg, M*MM); %#ok<MINV> %reorient lesion so midline is X=0        
%reslice to create a mirror image aligned in native space
P            = char([hdr.fname,',1'],[hdr_flip.fname,',1']);
flags.mask   = 0;
flags.mean   = 0;
flags.interp = 1;
flags.which  = 1;
flags.wrap   = [0 0 0];
flags.prefix = 'r';
spm_reslice(P,flags); 
delete(fname_flip); %remove flipped file
fname_flip = fullfile(pth,['rLR' nam ext]);%resliced flip file
%load lesion, blur 
hdrLesion = spm_vol(lesionImg); 
imgLesion = spm_read_vols(hdrLesion);
rdata = +(imgLesion > 0); %binarize raw lesion data, + converts logical to double
spm_smooth(rdata,imgLesion,4); %blur data
rdata = +(imgLesion > 0.1); %dilate: more than 20%
spm_smooth(rdata,imgLesion,8); %blur data
%now use lesion map to blend flipped and original image
hdr = spm_vol(anatImg); 
img = spm_read_vols(hdr);
hdr_flip = spm_vol(fname_flip); 
imgFlip = spm_read_vols(hdr_flip);
rdata = (img(:) .* (1.0-imgLesion(:)))+ (imgFlip(:) .* imgLesion(:));
rdata = reshape(rdata, size(img));
delete(fname_flip); %remove resliced flipped file
hdr_flip.fname = fullfile(pth,['e' nam ext]);%image with lesion filled with intact hemisphere
spm_write_vol(hdr_flip,rdata); 
intactImg = hdr_flip.fname;
%end entiamorphicSub()