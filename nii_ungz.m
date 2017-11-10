function Vout = nii_ungz (V, isCell, del)
%Unzips a .nii.gz image to .nii
%  V: image(s) to decompress
%  isCell: are results a cell strings or char strings?
%  del : (optional) if true delete original .gz image (FSL does not like
%         co-existing img.nii and img.nii.gz)
%Output: list of unzipped images
% Examples
%   nii_ungz('brain.nii.gz');
if ~exist('V','var') || isempty(V) %no files specified
 V = spm_select(inf,'^.*\.(gz|voi)$','Select gz files to decompress');
end;
if ischar(V), V = cellstr(V); end
Vout = {};
for i=1:numel(V)
  ref = deblank(V{i});
  [pth,nam,ext] = spm_fileparts(ref);
  if (length(ext)==3)  && min((ext=='.gz')==1)
    gunzip(ref);
    Vout ={Vout{:} fullfile(pth, [nam])};
    if exist('del','var') && del %del not specified
         delete(ref);
    end;
  elseif (length(ext)==4)  && min((ext=='.voi')==1)
    unz = gunzip(ref);
    [upth,unam,uext] = spm_fileparts(strvcat(unz)); %#ok<REMFF1>
    if isempty(uext) %if "file.voi" -> "file" then -> "file.nii"
        uext = '.nii';
        movefile(strvcat(unz),fullfile(upth, [unam uext])); %#ok<REMFF1>
    end;
    Vout ={Vout{:} fullfile(upth, [unam uext])};
  else
    Vout = {Vout{:} ref};
  end;
end; %for each file
if exist('isCell','var') && isCell, return; end;
Vout = strvcat(Vout); %#ok<REMFF1>
%end nii_ungz()