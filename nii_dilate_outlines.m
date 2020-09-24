function nii_dilate_outlines(fnms, voxs)
%Matryoshka doll effect
% fnms : file(s) to process
% voxs : distance of rings, e.g. [2, 4, 6] will create 2..4 and 4..6mm outlines
%Example
% nii_dilate_outlines('lesion.nii', [0 2])
if ~exist('fnms','var') %no files
    fnms = spm_select(inf,'image','Select images to dilate');
end
if ~exist('voxs','var') %no FWHM specified
    voxs = inputdlg({'Dilation in voxels (e.g. "2 4 6" for 2..4 and 4..6mm dilations)'},'Enter integers', [1],{'2 4 6'});
    voxs = str2num(voxs{1}); %#ok<NASGU>
end
if numel(voxs) < 2
    error('voxs must have at least 2 values (inner and outer boundary of outline)');
end
for i=1:size(fnms,1)
  fnm = deblank(fnms(i,:));
  [pth,nam,ext] = spm_fileparts(fnm);
  src = fullfile(pth,[nam ext]);
  hdr = spm_vol(src);
  if (voxs(1) == 0)
    inner = spm_read_vols(hdr);
  else
    inner = nii_dilate_vox(fnm, voxs(1));
  end
  for j = 2:numel(voxs)
    if (voxs(j) == 0)
        outer = spm_read_vols(hdr);
    else
        outer = nii_dilate_vox(fnm, voxs(j));
    end
    img = outer - inner;
    hdr.fname = fullfile(pth,[sprintf('d%g_%g',voxs(j-1), voxs(j)),  nam, ext]);
    hdr.dt(1)=2;%save binary data as bytes: uint8=2; int16=4; int32=8; float32=16; float64=64
    hdr.pinfo(1:2) = [1 0]; %scl_slope, scl_inter
    spm_write_vol(hdr, img);  
    inner = outer;
  end
end
%end nii_dilate_outlines()
