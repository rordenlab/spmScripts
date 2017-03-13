function fnm = nii_dti_clean(fnm, thresh);
%Detect and remove DTI volumes with substantial movement artifacts
%Warning: only for DTI: for fMRI we need to impute missing timepoints
% fnm : name of file to check and decimate if required
%Alternative
% fsl's eddy with outlier replacement (--repol)
%Examples
% nii_dti_clean; %use GUI
% nii_dti_clean('MUSC_2.7.nii')

if ~exist('fnm','var')
   %fnm = 'MUSC_2.7_AVE1_MCBI_40_dir_1007_V1_22.nii';
   fnm = spm_select(1,'^.*\.(gz|nii)$','Select DTI images to clean');
end
if ~exist('thresh','var')
   zThresh = 4;
end
[hdr,img, bvec,bval] = loadSub(fnm);
nvol = size(img,4);
nslice = hdr(1).dim(3);
if (nvol < 6) || (nslice < 5)
    error('%s requires DTI images with at least 7 volumes and 5 slices', mfilename);
end
%collapse X/Y dimension, so we have 3 dimensions (X*Y, Z (slices), volumes
img = reshape(img,size(img,1)*size(img,2),size(img,3),size(img,4));
stdev = zeros(nvol,1);
for v = 1: nvol
    dx = zeros(nslice-2,1);
    %mn = zeros(nslice-2,1);
    %mnTB = zeros(nslice-2,1);
    for s = 2: (nslice-1)
        mn = mean(img(:,s,v));
        mnTB = 0.5*(mean(img(:,s-1,v)+mean(img(:,s+1,v))));
        dx(s-1) = mn - mnTB;
    end;
    stdev(v) = std(dx)/mean(mean(img(:,:,v)));

end
stdev = stdev - mean(stdev(:)); %mean = 0
z = stdev ./ std(stdev(:));
%identify outliers
reject = zeros(nvol,1);
for v = 1: nvol
    if z(v) > zThresh
        fprintf('Deleting volume %d, z = %g\n', v, z(v));
        reject(v) = 1;
    elseif z(v) > (0.5 * zThresh)
        fprintf('Please visually inspect volume %d, z = %g\n', v, z(v));
    end
end
if sum(reject(:)) < 1
   fprintf('%s quitting: no volumes z>%g (max z=%g) in %s\n', mfilename, zThresh, max(z(:)), fnm);
   return;
end
fprintf('%s found %d volumes with z>%g in %s\n', mfilename, sum(reject(:)), zThresh, fnm);
%return raw image to 4D
img = reshape(img,hdr(1).dim(1), hdr(1).dim(2), hdr(1).dim(3), nvol);
%create decimated image
imgOK = img(:,:,:, ~reject(:));
nvolOK = sum(~reject(:));
%save image
hdrOK = hdr(1);
[p,n,x] = fileparts(hdr(1).fname);
n = [n '_' num2str(nvolOK)];
hdrOK.fname = fullfile(p, [n x]);
for v=1:nvolOK
    hdrOK.n(1)=v;
    spm_write_vol(hdrOK,imgOK(:, :, :, v));
end;
fnm = hdrOK.fname;
%save bvec/bval
bvecOK = bvec(~reject(:), :);
bvalOK = bval(~reject(:));
dlmwrite(fullfile(p,[ n, '.bval']), bvalOK','delimiter','\t');
dlmwrite(fullfile(p,[n, '.bvec']), bvecOK','delimiter','\t');
%end nii_dti_clean()


function [hd,im, bve,bva] = loadSub(fnm)
im = []; bve = []; bva = [];
[p,n,x] = spm_fileparts(fnm);
if (length(x)==3)  && min((x=='.gz')==1)
    fnm = char(gunzip(fnm));
    delnam = fnm;
    [p,n,x] = spm_fileparts(char(fnm));
else
    delnam = '';
end
hd = spm_vol(fnm); %input header
if hd(1).dt(1) == 128
    fprintf('Warning: skipping RGB image %s\n', fnm);
    if ~isempty(delnam), delete(delnam); end;
    return;
end
im = spm_read_vols(hd);%Input image
hd = hd(1);
if ~isempty(delnam), delete(delnam); end;
bve = loadTxt(fullfile(p,[n,'.bvec']), 3);
if isempty(bve), bve = zeros( size(im,4),3); end;
bva = loadTxt(fullfile(p,[n,'.bval']));
if isempty(bva), bva = zeros( size(im,4),1); end;
%end loadSub()

function txt = loadTxt(fnm, nRow)
txt = [];
%fprintf('%s\n', fnm);
if exist(fnm,'file') == 0, return; end;
fileID = fopen(fnm);
txt = cell2mat( textscan(fileID,'%f'));
fclose(fileID);
if exist('nRow', 'var')
    txt = reshape(txt, [numel(txt)/nRow nRow]);
end
%end loadTxt()