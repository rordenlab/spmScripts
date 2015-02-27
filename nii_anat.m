function nii_anat (Img, Outtext);
% Find normalized mm for MRIcron anat file
%  Img : nifti image before normalizeation
%  Outtext: name for output file [optional]
% Assumes MRIcron .anat and SPM _seg_inv_sn.mat files
%  e.g. 'img.nii' with 'img.anat' and 'img_seg_inv_sn.mat
% Depends on nii_map_coords.m
%Examples
%   nii_anat ('wAS_T1.nii');
%   nii_anat ('wAS_T1.nii', 'results.tab');

if nargin <1 %no image specified
 Img = spm_select(1,'image','Select images to convert');
end;
if nargin <2 %no text file
    Outtext = '';
end;

[pth,nam,ext] = spm_fileparts(Img);
anat = fullfile(pth,[nam '.anat']);
mat = fullfile(pth,[nam '_seg_inv_sn.mat']);

if ( (exist(Img) == 0) || (exist(anat) == 0) || (exist(mat) == 0) )
    fprintf('nii_anat unable to find required files %s %s %s\n',Img,anat,mat);
    return;
end;
[rowHeaders,vx_list] = readanatsub(anat);

if length(Outtext) > 0, myfile = fopen(Outtext ,'at'); end; 
rows = length(vx_list(:, 1));
if ( (length(rowHeaders(:)) ~= rows) || (length(vx_list(1, :)) ~= 3) || (rows < 1)   ) fprintf('Problem reading %s\n',anat); return; end;
for i = 1 : rows
	XYZ_mm =  vx_list(i, :)';
	[XYZ_mm XYZ_vx] = nii_map_coords(XYZ_mm, Img); % (XYZ_mm unaltered)
	[wXYZ_mm wXYZ_vx] = nii_map_coords(XYZ_vx, '', mat);
     %wXYZ_mm = XYZ_mm; % <- make a copy of original unwarped data....
    if length(Outtext) > 0
        fprintf(myfile,'%d\t%s\t%s\t%f\t%f\t%f\n',i,anat,strvcat(rowHeaders{i}),wXYZ_mm(1),wXYZ_mm(2),wXYZ_mm(3));
    else
        fprintf('%d\t%s\t%s\t%f\t%f\t%f\n',i,anat,strvcat(rowHeaders{i}),wXYZ_mm(1),wXYZ_mm(2),wXYZ_mm(3));
    end;
end;

if length(Outtext) > 0, fclose(myfile); end;



function [rowHeaders,num_list] = readanatsub(fileName)
%  Syntax to be used: [rowHeaders, num_list] =  readanatsub('filename.anat')
%  where 
%      rowHeaders will be a cell array
%      num_list will be a single-precision array
%  The tab-delimited input text-file should be formatted as follows:
%        F1  -13  16  -17
%        SY  3     4    5
%        VZ  4    -5   16
%  Adapted by Chris Rorden for mricron .anat files
%    based on Manu Raghavan,June 25, 2008
% [rowHeaders,num_list] = readanatsub('wAS_T1.anat')
%   rowHeaders{1} = 'F1'
%   num_list(1, :) = -13    16   -17
fid=fopen(fileName);
row = 1;
tline = fgetl(fid); % Get second row (first row of data)
while(1)
    tabLocs=findstr(char(9),tline); % find the tabs
    c = textscan(tline,'%s%f32%f32%f32%f32','Delimiter',char(9));
    rowHeaders{row} = c{1}; % Get column header
    for i=2:length(c)-1
        num_list(row,i-1) = c{i}; % Get numeric data
    end
    tline = fgetl(fid); % Go to next line in text file
    if(length(tline)==1)
        if(tline==-1) % Reached end of file, terminate
            break
        end
    else
        row = row+1;
    end        
end;
fclose(fid);