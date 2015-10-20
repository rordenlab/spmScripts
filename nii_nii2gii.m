function nii_nii2gii(filename, threshold, reduceFrac, smooth)
%convert NIfTI-format volume of voxels to a GIfTI-format mesh
% filename   : name of NIfTI image to convert
% threshold  : threshold for isosurface, e.g. "2" for all voxels brighter than 2
% reduceFrac : mesh simplification factor (0..1), e.g. 0.2 for 80% decimation of mesh
% smooth     : blur image
%Examples
% nii_nii2gii %use graphical interface
% nii_nii2gii('mx.nii', 4, 0.2, 1);
%Chris Rorden 2015, BSD license
% https://github.com/bonilhamusclab/MRIcroS/blob/master/%2BfileUtils/readVox.m

if ~exist('filename', 'var')  
   [A,Apth] = uigetfile({'*.nii;*.hdr;*.gz;'},'Select NIFTI image', 'MultiSelect', 'off');
   filename = strcat(Apth,char(A));
end;
if ~exist(filename,'file'), error('Unable to find image'); end;
if ~exist('smooth', 'var')    
    prompt = {'Enter threshold:','Enter reduction factor','Smooth amount (0=none, 1=1voxel, 2=2voxel)'};
    dlg_title = 'Meshify options';
    num_lines = 1;
    def = {'2','0.3','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    threshold = str2double(answer{1});
    reduceFrac = str2double(answer{2});
    smooth = str2double(answer{3});
end;
%decompress image if saved as nii.gz
[p,n,x] = fileparts(filename);
if strcmpi(x,'.gz')
    filename = char(gunzip(filename));
end
%read image
hdr = spm_vol(filename);
img = spm_read_vols(hdr);
if strcmpi(x,'.gz')
    delete(filename ); %FSL hates filename.nii co-existing with filename.nii.gz
end
%apply threshold
img(isnan(img)) = 0; %remove not-a-numbers
if smooth > 0 %blur image prior to edge extraction
    smoothD = round(smooth)*2+1; %smooth MUST be an integer
    fprintf('Applying gaussian smooth with %d voxel diameter\n',round(smoothD));
    img = smooth3(img,'gaussian',round(smoothD), round(smoothD) * 0.2167);
end;
if (threshold < 0)
    img = -img;
    threshold = -threshold;
end
if (threshold > max(img(:))), error('No voxels survive threshold, max: %g', max(img(:))); end;
FV = isosurface(img,threshold);
if reduceFrac < 1
    FV = reducepatch(FV,reduceFrac);
end;
%next: isosurface swaps the X and Y dimensions! size(Vol)
FV.vertices =  FV.vertices(:,[1:0,2,2:1,1,3:end]);
%convert from voxels to mm
vx = [ FV.vertices ones(size(FV.vertices,1),1)];
vx = mtimes(hdr.mat,vx')';
FV.vertices = vx(:,1:3);
%something seems to reverse faces!
FV.faces = fliplr(FV.faces); 
%save data
outname = fullfile(p, [n, '.gii']);
save(gifti(FV),outname,'GZipBase64Binary');
%nii_nii2gii()
