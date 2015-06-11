function nii_rescue_par (parname)
%Fix broken PAR/REC files where scan was interrupted mid-volume
% dtiNii: name of bvec file(s), e.g. img.bvec
%Examples
% nii_rescue_par
% nii_rescue_par('fMRI.PAR');
%Chris Rorden 11 June 2015, BSD License


if ~exist('parname','var') %file not specified
   [A,Apth] = uigetfile({'*.PAR';'*.*'},'Select PAR file(s)', 'MultiSelect', 'off');
   parname = strcat(Apth,char(A));
end;
[pth,nam] = fileparts(parname);
recname = fullfile(pth, [nam, '.REC']);
if ~exist(recname,'file'), error('Unable to find file %s', recname); end;
%step 1 read par header
parlines = textread(parname, '%s', 'delimiter', '\n', 'whitespace', '');
mxSlice = 0;
mxVol = 0;
numSlice = 0;
for i = 1: numel(parlines)
    str = char(parlines(i));
    if isempty(str) || (str(1) == '.') || (str(1) == '#'), continue; end; 
    vals = strread(str,'%f');
    sl = vals(1);
    vol = vals(3);
    if sl > mxSlice, mxSlice = sl; end;
    if vol > mxVol, mxVol = vol; end;
    numSlice = numSlice + 1;
end
%check that we can restore this image
if numSlice == (mxVol * mxSlice), error('Volumes divisible by slices: no need to fix'); end;
if ((mxVol-1) * mxSlice) > numSlice, error('PAR file borked in unknown manner'); end;
recDir = dir(recname);
if rem(recDir.bytes,numSlice) ~= 0, error('Unable to determine slice size'); end;
%write new file...
sliceBytes = recDir.bytes/numSlice;
parnameOut = fullfile(pth, ['x', nam, '.PAR']);
recnameOut = fullfile(pth, ['x', nam, '.REC']);
parOutFID = fopen(parnameOut,'w');
recFID = fopen(recname, 'r');
recOutFID = fopen(recnameOut, 'wb','native');
for i = 1: numel(parlines)
    str = char(parlines(i));
    if isempty(str) || (str(1) == '.') || (str(1) == '#')
        fprintf(parOutFID,'%s\n',str);
        continue; 
    end; 
    vals = strread(str,'%f');
    A = fread(recFID, sliceBytes, 'uint8');
    if vals(3) >= mxVol, 
        continue; 
    end;
    fwrite(recOutFID,A, 'uint8');
    fprintf(parOutFID,'%s\n',str);
end
fclose(parOutFID);
fclose(recFID);
fclose(recOutFID);