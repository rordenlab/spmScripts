function nii_newseg (P,norm, Template,T2);
%SPM8 new segment with Chris Rordens extended TPM and tissue cleanup
%  P: Images to normalize
%  norm : [optional] Will results be saved in normalized (true) or native (false) space?
%  template: [optional] name of template to use (defaults to extended FOV eTPM.nii)
%
%Example
%  nii_newseg('T1.nii');

cleanup = 1;

if nargin <1 %no files
 P = spm_select(inf,'image','Select images for new segment');
end;
if nargin <2 %norm not specified... do not normalize data
    norm = false;
end;
if nargin <3 %no Template
 %Template = fullfile(spm('Dir'),'toolbox','Seg','TPM.nii');%SPM8 default template
 Template = fullfile(fileparts(which(mfilename)),'eTPM.nii');
 %Template = spm_select(1,'image','Select template for new segment');
end;

if norm
    n = [1 0];
    w = [1 0];
else
    n = [1 1];
    w = [0 0];
end;

spm('defaults','fmri');
spm_jobman('initcfg');
[ptht,namt,extt] = spm_fileparts(deblank(Template(1,:)));
Tem = [ptht,filesep,namt,extt];
spm_jobman('initcfg');

for i=1:size(P,1)
  ref = deblank(P(i,:));
  [pth,nam,ext] = spm_fileparts(ref);
  matlabbatch{1}.spm.tools.preproc8.channel.vols = {ref};
  matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
  matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
  matlabbatch{1}.spm.tools.preproc8.channel.write = [0 1];
  if nargin>3 %T2 specified
    ref2 = deblank(T2(i,:));
    [pth2,nam2,ext2] = spm_fileparts(ref2);
	fprintf('Using %s to segment T1 and T2 images: %s %s\n', Template, ref, ref2);
    matlabbatch{1}.spm.tools.preproc8.channel(2).vols = {ref2};
	matlabbatch{1}.spm.tools.preproc8.channel(2).biasreg = 0.0001;
	matlabbatch{1}.spm.tools.preproc8.channel(2).biasfwhm = 60;
	matlabbatch{1}.spm.tools.preproc8.channel(2).write = [0 0];
  else
      	fprintf('Using %s to segment T1 %s\n', Template, ref);
  end;
  matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[Tem,',1']};
  matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
  matlabbatch{1}.spm.tools.preproc8.tissue(1).native = n;
  matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = w;
  matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[Tem,',2']};
  matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
  matlabbatch{1}.spm.tools.preproc8.tissue(2).native = n;
  matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = w;
  matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[Tem,',3']};
  matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
  matlabbatch{1}.spm.tools.preproc8.tissue(3).native = n;
  matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = w;
  matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[Tem,',4']};
  matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
  matlabbatch{1}.spm.tools.preproc8.tissue(4).native = n;
  matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = w;
  matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[Tem,',5']};
  matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
  matlabbatch{1}.spm.tools.preproc8.tissue(5).native = n;
  matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = w;
  matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[Tem, ',6']};
  matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
  matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
  matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
  matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
  matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
  matlabbatch{1}.spm.tools.preproc8.warp.write = w;
  spm_jobman('run',matlabbatch);
  if cleanup 
  	nii_cleanup5([pth,filesep,'c1',nam,ext],[pth,filesep,'c2',nam,ext],[pth,filesep,'c3',nam,ext], [pth,filesep,'c4',nam,ext],[pth,filesep,'c5',nam,ext],[pth,filesep,'m',nam,ext]); 
  end;
end;%for each image