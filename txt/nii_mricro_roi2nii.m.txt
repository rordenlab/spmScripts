function nii_mricro_roi2nii (roinames, imgnames)
%Convert old MRIcro ROI format to NIfTI format
% roinames: name(s) of .roi file
% imgnames: name(s) of NIfTI format image roi was drawn on

if ~exist('roinames','var') %no files specified
    roinames = spm_select(inf,^.*\.roi$,'Select Region(s) of Interest');
end
if ~exist('imgnames','var') %no files specified
 imgnames = spm_select(inf,'image','Select background images (same order)');
end
if (size(imgnames,1) ~= size(roinames,1))
    error('You must select one image for each ROI');
end
for i=1:size(imgnames,1)
    imgname = deblank(imgnames(i,:));
    roiname = deblank(roinames(i,:));
    hdr = spm_vol(imgname); %load header 
    img = spm_read_vols(hdr); %load image data
    oimg = zeros(size(img));
    %load raw binary data
    fid=fopen(roiname,'rb'); % 
    raw = fread(fid, inf, 'uint32');
    fclose(fid);
    %read data
    sliceVox = hdr.dim(1)* hdr.dim(2);
    max12bit = 2^12-1;
    max15bit = 2^15-1;
    max16bit = 2^16-1;
    if sliceVox > max16bit
        %12 bit runs and 20 bit offsets
        signature = 1;
        maxRun = max12bit;
    else
        signature = 0;
        maxRun = max16bit;
        %16 bit runs and 16 bit offsets  
    end
    if (bitand(bitshift(raw(1),-15),1) ~= signature) 
        %each roi file starts with a signature bit:
        %  0 indicates the ROI slice size can be stored with 16 bits (<65535 voxels)
        %  1 indicates larger than 16 bit slice size
        error('The roi and image do not match (1: signature incorrect)');
    end
    nraw = numel(raw);
    pos = 1;
    while true
        runs = bitshift(raw(pos),-17)-1;
        slice = bitand(raw(pos),max15bit);
        %fprintf('slice = %d runs = %d  %d\n',slice, runs);
        if slice > hdr.dim(3)
            error('The roi and image do not match (2: slice out of range)'); 
        end
        slice = (slice-1) * sliceVox;
        pos = pos + 1;
        for r = 1:runs
            if pos > nraw, break; end;
            len = bitand(bitshift(raw(pos),-16),maxRun);
            offs = bitand(raw(pos),max16bit);
            if signature == 1;
                offs = offs + bitshift(bitshift(raw(pos),-28),+16);
            end
            if (offs+len) > (sliceVox+1)
                error('The roi and image do not match (3: roi exceeds slice size)'); 
            end
            %fprintf(' len = %d offs = %d\n',len,offs);
            oimg(slice+offs:slice+offs+len-1) = 1;
            pos = pos + 1;
        end %for each run
        if pos > nraw, break; end;
    end
    if (sum(hdr.mat(1:3)) < 0)
        oimg = flipdim(oimg,1);
    end
    %[pth, nm, ~] = fileparts(roiname);
    [pth, nm, ~] = fileparts(imgname);
    hdr.fname = fullfile(pth, ['l' nm '.nii']);
    spm_write_vol(hdr,oimg);
end