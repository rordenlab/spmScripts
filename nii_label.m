function nii_label (fnm, threshold)
%Threshold volume, then label each cluster in the volume
% fnm  : name of image
% threshold : intensity threshold, e.g. if 2.0 clusters where intensity exceeds 2.0 will survive.
%Examples
% nii_label %prompt user for files and values
% nii_label('img.nii',2.5);
%A wrapper for Jesper Andersson's elegant spm_bwlabel

%provide user interface if parameters not specified
if ~exist('fnm','var')
 fnm = spm_select(1,'image','Select images to threshold');
end;
if ~exist('threshold','var')
    answer = inputdlg({'Intensity threshold (e.g. 2 for clusters brighter than 2)'}, 'Set threshold', 1,{'2'});
    threshold = str2double (cell2mat(answer(1)));
end
%load image
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
%binarize data
if threshold < 0 
    img(img > threshold) = 0;
else
    img(img < threshold) = 0;
end
img(img ~=0) = 1;
if sum(img(:)~= 0) == 0
    error('No voxels survive this threshold');
end
%report results
[imgLabel,nLabel] = spm_bwlabel(img,18);
fprintf('%s has %d voxels in %d clusters that exceed threshold %f\n', hdr.fname,sum(img(:)~= 0),nLabel, threshold );
%report coordinates, save txt file
[pth nm ext] = spm_fileparts(fnm);
fid = fopen( fullfile(pth, ['roi_' nm '.txt']), 'wt' );
for label = 1: nLabel 
    imgL = imgLabel;
    imgL(imgL ~= label) = 0;
    vx = centerOfMassSub (imgL);
    XYZ_mm = hdr.mat * [vx 1]';
    coord = XYZ_mm(1:3)';
    txt=mat2str((coord));
    fprintf('Label %d has %d voxels with a center of mass %s\n',label,sum(imgL(:)~= 0), txt);
    txt=mat2str(round(coord));
    txt(txt==' ')='_';
    txt(txt=='[')='';
    txt(txt==']')='';
    fprintf( fid, '%d|%s\n', txt);
end
fclose(fid);
%save image data
hdr.fname = fullfile(pth, ['roi_' nm ext]);
spm_write_vol(hdr,imgLabel);
%end nii_label()

function coord = centerOfMassSub (A)
% Jared Wells Center of Mass, BSD licence http://www.mathworks.com/matlabcentral/fileexchange/41675-center-of-mass
sz = size(A);
nd = ndims(A);
M = sum(A(:));
coord = zeros(1,nd);
if M==0
    coord = [];
else
    for ii = 1:nd
        shp = ones(1,nd);
        shp(ii) = sz(ii);
        rep = sz;
        rep(ii) = 1;
        ind = repmat(reshape(1:sz(ii),shp),rep);
        coord(ii) = sum(ind(:).*A(:))./M;
    end
end
%end centerOfMassSub()
