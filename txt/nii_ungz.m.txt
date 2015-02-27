function Vout = nii_ungz (V)
% Unzips a .nii.gz image to .nii
% Input: list of images
% Output: list of unzipped images
%
% Example
%   nii_ungz('brain.nii.gz');

if nargin <1 %no files
 V = spm_select(inf,'^.*\.gz$','Select gz files to decompress');
end;

Vout = {};
for i=1:size(V,1)
  ref = deblank(V(i,:));
  [pth,nam,ext] = spm_fileparts(ref);
  if (length(ext)==3)  && min((ext=='.gz')==1) 
    gunzip(ref);
    Vout ={Vout{:} fullfile(pth, [nam])};
    
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
Vout = strvcat(Vout); %#ok<REMFF1>
