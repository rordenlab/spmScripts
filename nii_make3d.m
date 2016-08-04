function nii_make3d
%simple script to generate 4d nifti image - useful for testing other routines
% http://paulbourke.net/geometry/borg/

inname = fullfile(spm('Dir'),'apriori','brainmask.nii');
if ~exist(inname, 'file')
    inname =  fullfile(spm('Dir'),'canonical','avg152T1.nii');
    if ~exist(inname, 'file')
        fprintf('%s error: unable to find image named %s\n', mfilename,inname);
        return;
    end;
end
hdr = spm_vol([inname,',1']); 
img = spm_read_vols(hdr);
nam = 'test3d.nii';
[nX nY nZ] = size(img);
img = zeros(nX,nY,nZ); 
freq =7;
for z = 2 : (nZ-1) 
	z1 = freq*z/nZ;
	for y = 2 : (nY-1) 
		y1 = freq*y/nY;
		for x = 2 : (nX-1) 
			x1 = freq*x/nX;			
			img(x,y,z)=pi+sin(x1*y1)+sin(y1*z1)+sin(z1*x1);
		end;
	end;
end;
hdr.fname = nam;  
hdr.pinfo = [0;0;0];
spm_write_vol(hdr,img);
