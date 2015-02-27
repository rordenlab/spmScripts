function nii_make4d
%simple script to generate 4d nifti image - useful for testing other routines
nvol = 64;
nam = 'test4dx.nii';
dt = 2; %desired format of converted image uint8=2; int16=4; int32=8; float32=16; float64=64
bigEndian = false; %big or little endian format
tr= 1; %e.g. 2.2 if 2.2s per volume
drift = 1; %e.g. if 0.1 then final volume 0.1 brighter than first
drift = drift/nvol;
%we traditionally try to preserve 0.01-0.1Hz signals 
slowestFreq = 0.00390625;
highestFreq = 0.25;
inname = fullfile(spm('Dir'),'apriori','brainmask.nii');
if ~exist(inname, 'file')
    fprintf('%s error: unable to find image named %s\n', mfilename,inname);
    return;
end
hdr = spm_vol([inname,',1']); 
img = spm_read_vols(hdr);
img(img(:) == 0) = nan; 
nZ = size(img,3);
%create one sine wave for each slice
fprintf('volume TR is %f\n',tr);
t=0:tr:(tr*nvol); %onset time of each volume
%t = t + 1; %all frequencies are zero at time zero
y=zeros(nZ,length(t));
freqInc = (highestFreq-slowestFreq) / (nZ-1);
for z=1:nZ
        f = slowestFreq + ((z-1)*freqInc);
        y(z,:)=sin(2*pi*f*t);
        fprintf('slice %d has a frequency of %f\n',z,f);
end;
imgVol1 = img(:, :, :, 1); %first volume of data
%embed column in volume so all slices have clear signal
for zi=1:nZ
    for yi=10:20
        for xi=10:20
            imgVol1(xi,yi,zi) = 1;
        end
    end;
end;
img = zeros(size(imgVol1,1),size(imgVol1,2),size(imgVol1,3),nvol);
for vol=1:nvol
    imgMod = imgVol1;
    for z=1:nZ
        imgMod(:,:,z) = imgVol1(:,:,z)*y(z,vol);
    end;
    imgMod = imgMod + imgVol1; %shift baseline away from zero
    imgMod = imgMod + drift*vol; %add drift
    linearTrend = 1+ drift*vol; %make a block that only shows linear trend
    for yi=25:35
        for xi=10:20
            imgMod(xi,yi,:) = linearTrend;
        end
    end;
    for yi=45:55 %constant band
        for xi=10:20
            imgMod(xi,yi,:) = 1;
        end
    end;
    img(:,:,:,vol) = imgMod(:,:,:);
end;
%normalize intensity
img = img - min(img(:)); %translate so minimum is zero
img = img/max(img(:)); %scale so maximum is one
%set nans (regions outside brain) to zero
img(isnan(img)) = 0; 
%set output
pth = pwd;
hdr.fname   = fullfile(pth,nam );
if exist(hdr.fname, 'file')==2
  delete(hdr.fname);
end
hdr.private.timing.toffset= 0;
hdr.private.timing.tspace= tr;
%next: save image with desired precision
%hdr.dt    =[dt,0];
hdr.dt    =[dt,bigEndian]; %hdr.dt(2) determines little (hdr.dt(2)=0) or big (hdr.dt(2)=1) endian
if isnan(spm_type(hdr.dt,'intt'))
    error('Unknown datatype %d\n',dt);
end
if (spm_type(hdr.dt,'intt')) %integer input
    %from spm_file_merge
    fprintf('Saving %s as integers with %d-bits per voxel\n',hdr.fname, spm_type(hdr.dt,'bits'));
    dmx  = spm_type(hdr.dt,'maxval');
    dmn  = spm_type(hdr.dt,'minval');
    mx = max(img(:));
    mn = min(img(:));
    if dmn < 0
    	sf = max(mx/dmx,-mn/dmn);
    else
    	sf = mx/dmx;
    end
    hdr.pinfo = [sf;0;0];
else
    fprintf('Saving %s as real-numbers (floating point) with %d-bits per voxel\n',hdr.fname, spm_type(hdr.dt,'bits'));
    hdr.pinfo = [1;0;0]; 
end
fprintf('Ignore potential warnings that data types do not match\n');%this may complain is dt does not match hdr.private.dat.dtype = [];
for vol=1:nvol
    hdr.n(1)=vol;
    
    spm_write_vol(hdr,img(:, :, :, vol));
end;