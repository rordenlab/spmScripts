We assessed the most reliable links based on the intact (right) hemisphere of all 145 participants. Specifically, for each individual we computed the mean connectivity of each of the 53 right hemisphere gray matter regions to all the other 52 regions. The mean and variability of these measures allowed us to compute a Z-score, where a high z-score reflects a relatively large mean effect to a relatively small variability in mean. We then rank ordered these Z-scores to determine the most reliable connections using our methods.

  Here is the mean data for 145 stroke patients. To view the data with Matlab type:
l = load('fiberQA.mat')
            z: [53x53 double]
       labels: [53x117 char]
      numSubj: 145
    labelUsed: [189x1 logical]
The values are as follows:
 z: this is a matrix of z-scores for the 53 regions connected to all the other 53 regions
 labels: names of the 53 regions (all right hemisphere gray matter)
 numSubj: number of individuals who contributed
 labelUsed: this vector maps the 53 gray matter regions to the complete set of 189 JHU regions


Here is a sample script to show you how to find the strongest link:

l = load('fiberQA.mat');
[mx, indx] = max(l.z(:));
[mxRow,mxCol] = ind2sub(size(l.z),indx);
fprintf('strongest link is "%s"*"%s"\n',deblank(l.labels(mxRow,:)), deblank(l.labels(mxCol,:)))

Which reports

strongest link is "92|RedNc_R|red nucleus right|1"*"100|Midbrain_R|midbrain right|1"
