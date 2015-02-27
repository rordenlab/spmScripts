function nii_cleanup5_batch;
%clean segmentation partitions for 5 tissue types, optionally create image for volume rendering

T1 = spm_select(inf,'image','Select T1 images that were used for normalization');

for i=1:size(T1,1), %repeat for each image the user selected
    [p,n,x] = spm_fileparts(deblank(T1(i,:)));
    fprintf('%d of %d\n',i,size(T1,1))
    nii_cleanup5(fullfile(p,['rc1' n x]),fullfile(p,['rc2' n x]),fullfile(p,['rc3' n x]),fullfile(p,['rc4' n x]),fullfile(p,['rc5' n x]) ); 
end;

