function nii_norm_linear (V);
%SPM's "Normalise" function with only linear transforms
% By default, images are resliced to 1mm isotropic
%
% Example
%   nii_norm_linear('C:\dir\img.nii');

if nargin <1 %no files
 V = spm_select(inf,'image','Select images to coreg');
end;
[pth,nam,ext] = spm_fileparts(deblank(V(1,:)));
src = fullfile(pth,[nam ext]);
if (exist(src) ~= 2) 
 	fprintf('nii_norm_linear error: unable to find template image %s.\n',src);
	return;  
end;
%ref = fullfile(spm('Dir'),'templates','MNI152lin_T1_1mmlr.nii');
ref = fullfile(spm('Dir'),'templates','MNI152lin_T1_1mmlr.nii');
if (exist(ref) ~= 2) 
 	fprintf('Error: unable to find template image %s.\n',ref);
	return;  
end;

matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {[src ,',1']};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {[src ,',1']};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = {[ref,',1']};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = Inf;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
%matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [  -90 -126  -72;  90   90  108];
%matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [2 2 2];

matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [  -90.5 -126.5  -72.5;  90.5   90.5  108.5];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);