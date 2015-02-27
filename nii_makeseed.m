function nii_makeseed (V, Radius,Mask);
% Finds peak for each image V, generates image with sphere at peak
% V: Image[s] to create seed maps
% PeakRadius: Voxels radius for seed size at peak
% Mask: (optional). list of mask image[s] - constrains peak search to mask
%
%Example
%  nii_makeseed('brain.nii');
%  nii_makeseed('brain.nii.gz',4);
%  nii_makeseed('rain.nii.gz',4,'post.nii.gz');
%  nii_makeseed('rain.nii',4,strvcat('post.nii','pre.nii'));
 
if nargin <1 %no files
 %V = spm_select(inf,'image','Select images to average');
 V = spm_select(inf,'^.*\.(gz|voi|img|nii)$','Select image(s) to make peak seed(s)');
end
if nargin <2 %Radius not specified
 Radius = 8; %8 voxels
end;
for i=1:size(V,1)
    %load the image
	ref = deblank(V(i,:));
    ref = nii_ungz(ref);
	[pth,nam,ext] = spm_fileparts(ref); 
	fnm = fullfile(pth,[nam ext]);
	VO = spm_vol(fnm);
 	img = spm_read_vols(VO);
    if nargin >2 %Mask the image
        for j=1:size(Mask,1)
            %load mask
            mname = deblank(Mask(j,:));
            mname = nii_ungz(mname);
            [mpth,mnam,mext] = spm_fileparts(mname); 
            Vm = spm_vol(mname);
            m = spm_read_vols(Vm);
            mn = min(m(:));
            %make masked version of image
            mimg = img; 
            mimg((m==mn)) = 0; 
            VO.fname = fullfile(pth, ['a' mnam  ext]);
            findpeak(mimg,Radius,VO);
        end; %for each mask: j
    else %else no mask specified
        VO.fname = fullfile(pth, ['a' nam  ext]);
        findpeak(img,Radius, VO); 
    end; %if..else no mask specified
end;

function findpeak (img,  Radius, VO);
    %find the peak intensity
	[peakht peakloc] = max(img(:));
 	[p(1) p(2) p(3)] = ind2sub(size(img), peakloc);
    numpeak = length(find(img==peakht));
    if numpeak > 1
        fprintf('Warning: there are %d voxels with the peak intensity %f\n',numpeak, peakht);
    else
          fprintf('Peak intensity is voxel at row*column*slice %d*%d*%d = %f\n',p(1), p(2), p(3), peakht);      
    end;
    fprintf(' Sphere with radius of %f\n',Radius); 
	%draw a sphere centered at the peak
    d = size(img);
 	lo = p-Radius;
    hi = p+Radius;
    lo((lo<1)) = 1;
    if hi(1)>d(1) hi(1)=d(1); end;
    if hi(2)>d(2) hi(2)=d(2); end;
    if hi(3)>d(3) hi(3)=d(3); end;
	img(:) = 0;    
    	for z = lo(3) : hi(3) 
        	for y = lo(2) : hi(2) 
            		for x = lo(1) : hi(1)
                		%compute distance from origin using Pythagorean theorem
				dx = sqrt((p(1)-x)^2+(p(2)-y)^2+(p(3)-z)^2 );
                		if (dx<Radius) img(x,y,z)=1; end;
            		end;%X
        	end;%Y
	end;%Z
    %save the image
	img(isnan(img)) = 0;
	spm_write_vol(VO,img);