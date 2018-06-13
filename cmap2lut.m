function cmap2lut (fnm)
%convert fsleyes text format cmap color maps to ImageJ/MRIcron lut binary format
% fnm : file to convert
%Example
% cmap2lut('copper.cmap');

fnm = 'copper.cmap'
if ~exist('fnm','var')  
   [files,pth] = uigetfile({'*.cmap;';'*.*'},'Select the cmap[s]', 'MultiSelect', 'on'); 
   fnms = cellstr(files); 
else
    [pth,nam, ext] = fileparts(fnm);
    fnms = cellstr([nam, ext]); 
end;

for i=1:size(fnms,2) %apply to image(s)
    fnm = fullfile(pth,deblank(fnms{i}));
    %read text
    fid = fopen(fnm,'r');
    RGB = fscanf(fid,'%f %f %f');
    fclose(fid);
    RGB = reshape(RGB,3,numel(RGB)/3)';
    %we need precisely 256 entries
    i = 1;
    RGB256 = [];
    RGB6 = [];
    for i = 1 : 3
        RGB256 = [RGB256, interp1(linspace(0,1,size(RGB,1)), RGB(:,i)', linspace(0,1,256))'];
        RGB6 = [RGB6, interp1(linspace(0,1,size(RGB,1)), RGB(:,i)', linspace(0,1,6))'];

    end
    RGB6 = round(RGB6 * 255)
    return;
    %save binary
    [p,n] = fileparts(fnm);
    fid = fopen(fullfile(p,[n,'.lut']),'wb');
    fwrite(fid,RGB256*255,'uchar');
    fclose(fid);
end
%end cmap2lut()






