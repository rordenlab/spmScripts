function nii_reorient(fnms, M);
%Rotate NIfTI images 
% fnms : names of files
% m : matrix;
%Examples
% M = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %LR mirror
% nii_reorient('LR.nii', M);

fprintf('Warning %s reorients images by changing the NIfFTI sform\n', mfilename);
if ~exist('fnms','var') || isempty(fnms)
    %fnms = 'avg152T1.nii';
	fnms = spm_select(inf,'image','Select image[s] for NaN removal'); 
end
if ~exist('M','var')
	%M = [1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1]; %mirror
   prompt={'Enter rotation matrix'};
   name='Rotations';
   numlines=1;
   %defaultanswer={'-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1'};
   defaultanswer={'-1 0 0; 0 1 0; 0 0 1'};
   answer=inputdlg(prompt,name,numlines,defaultanswer);
   M = str2num(answer{1})
end
if (size(M,1) == 3) && (numel(M)==9) 
   %user provided a 3x3 matrix, pad to 4x4
   t = M;
   M = eye(4);
   M(1:3,1:3) = t;
end
if (size(M,1) ~= 4) && (numel(M)~=16) 
   error('Rotation matrix M must be a 3x3 or 4x4 matrix'); 
end
%process files
for i=1:size(fnms,1)
    fnm = deblank(fnms(i,:));
    hdr = spm_vol(fnm);
    img = spm_read_vols(hdr);
    rM = M*spm_get_space(fnm);
    spm_get_space(fnm,rM);
    [pth nm ext] = spm_fileparts(fnm);
    matnam = fullfile(pth, [nm, '.mat']);
    if exist(matnam, 'file')
        delete(matnam); %single transform for all volumes
    end;
end;
