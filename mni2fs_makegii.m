function mni2fs_makegii (mnivol, hem, smoothdata, InterpMethod)
%Converts NIfTI image in MNI space to GIfTI image in freesurfer space
% mnivol = filename of NIfTI file to convert, e.g. 'c:\AudMean.nii'
% hem = hemisphere to plot, either left ('lh') or right ('rh')
% smoothdata = amount to blur data, 0 for none
% InterpMethod = can be either 'nearest' 'linear' 'cubic' 'spline'
%Examples
% mni2fs_makegii('./examples/AudMean.nii','lh',0, 'linear');
% mni2fs_makegii('./examples/HOA_heschlsL.nii','lh',0, 'linear')
% mni2fs_makegii('./examples/motor_4t95.nii.gz','lh',0, 'linear')
% mni2fs_makegii; %use graphical interface


if ~exist('mnivol','var') %filename not specified: graphical interface
    [A,Apth] = uigetfile({'*.nii;*.hdr;*.gz;'},'Select NIFTI image in MNI space', 'MultiSelect', 'off');
    mnivol = strcat(Apth,char(A));
end
if ~exist('InterpMethod','var') %some/all arguments not specified: graphical interface
    prompt = {'Enter hemisphere ("lh" or "rh"):','Smoothing (0 = none):','InterpMethod ("nearest" "linear" "cubic" "spline")'};
    dlg_title = 'Values for conversion';
    num_lines = 1;
    def = {'lh','0','linear'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    hem = char(answer{1});
    if (numel(hem) < 1) || strcmpi(hem(1),'r')
        hem = 'rh'; 
    else
        hem = 'lh';
    end;
    smoothdata = str2double(answer{2});
    InterpMethod = char(answer{3});
end;

thisfolder = fileparts(mfilename('fullpath'));
addpath(genpath(thisfolder)) % will add all subfolders and dependencies
if ischar(mnivol)
    [p,n,x] = fileparts(mnivol);
    if strcmpi(x,'.gz') %decompress FSL-style nii.gz
        mnivol = char(gunzip(mnivol));
        NII = load_untouch_nii(mnivol);
        delete(mnivol); %remove decompressed file
    else
        NII = load_untouch_nii(mnivol);
    end
    gii_fn = fullfile(p,[n,'.gii']);
elseif isstruct(mnivol)
    NII = mnivol;
    gii_fn = 'vertexColors.gii';
end
if isinteger(NII.img) % Convert NII image to single
    NII.img = single(NII.img);
end
%NII.img(isnan(NII.img)) = 0;
if smoothdata > 0
    disp('Smoothing Volume')
    NII.img = smooth3(NII.img,'gaussian',15,smoothdata);
end

%Load all vertices
surfrender_fn = fullfile(thisfolder,['/surf/' hem '.pial.surf.gii']);
g = gifti(surfrender_fn);
V = g.vertices;
%%%
Tmni = [NII.hdr.hist.srow_x; NII.hdr.hist.srow_y; NII.hdr.hist.srow_z; 0 0 0 1];
sz = size(NII.img);
[X, Y, Z] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
X = X(:)-1; Y = Y(:)-1; Z = Z(:)-1;
L = length(X);
XYZmni = [X Y Z];
XYZmni = [XYZmni ones(L,1)]*Tmni'; % transform from voxel to mni space
X = reshape(XYZmni(:,1),sz);
Y = reshape(XYZmni(:,2),sz);
Z = reshape(XYZmni(:,3),sz);
% Get vertices in MNI space
load(fullfile(thisfolder,'surf/transmats.mat'),'Tfstovox_rcor','Trsvoxtomni_rcor');
V = [V ones(length(V),1)]*Tfstovox_rcor'*Trsvoxtomni_rcor';
Vsurf = interpn(X,Y,Z,NII.img,V(:,1),V(:,2),V(:,3),InterpMethod);
min(Vsurf(:))
max(Vsurf(:))

if (min(Vsurf(:)) == max(Vsurf(:)) ) 
    fprintf('No variability in output intensity (all voxels = %g). Input varies from %g..%g\n', min(Vsurf(:)), min(NII.img(:)), max(NII.img(:)) ); 
    return;


end;
% if true % Save check to file
%  Vmnivox = round(V*inv(Tmni)');
%  for ii = 1:length(Vmnivox)
%      NII.img(Vmnivox(ii,1),Vmnivox(ii,2),Vmnivox(ii,3)) = NII.img(Vmnivox(ii,1),Vmnivox(ii,2),Vmnivox(ii,3))+0.2;
%  end
%  save_untouch_nii(NII,'test_surf.nii');
% end
% save the data
g = gifti;
g = subsasgn(g, struct('type','.','subs','cdata'), Vsurf);
save(g,gii_fn,'GZipBase64Binary');