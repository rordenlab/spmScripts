function t2_contrast;
%plots T2-decay curves for two tissues as well as their relative contrast
%http://en.wikipedia.org/wiki/Relaxation_(NMR)
%http://en.wikipedia.org/wiki/Spin-spin_relaxation_time

%Note T2 depends on tissue and field strength:
% Table 2.1 of Magnetic resonance imaging of the brain and spine, Volume 1, By Scott W. Atlas page 27
%T2 (ms)
%         Gray     White	CSF     Blood(Venous)   Blood(Arterial)   Fat
%   1.5T  100      80       2200    150             235                165 
%   3.0T  85       77       2200    80              165                133
% Table 11.1 of Functional Neuroimaging in Exercise and Sport Sciences, Boecker et al.,  page 243
%T2 (ms)
%         cGray    White	CSF Blood(a)
%   1.5T  80       80       500  256
%   3.0T  70       60       500  121
%   7.0T  60       40        60   50
% for blood see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3752414/
%  both venous and arterial blood T2 decrease exponentially with B0 

x = 1:1:300;
xyGM = exp(-x./85);
xyFat = exp(-x./133);
xyContrast = abs(xyGM-xyFat);
p =plot(x,xyGM,x,xyFat,x,xyContrast);
xlabel('TE (ms)');
ylabel('Signal');
legend('Gray Matter','Fat','Relative Contrast');

set(p,'LineWidth',2);
set(gcf,'Color',[1 1 1]);

%next: report TE that has optimal contrast
[mx, i] = max(xyContrast);
title( ['Optimal GM/Fat T2 contrast at ',int2str(x(i)), 'ms (3.0T)']);