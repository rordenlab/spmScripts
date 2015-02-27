function t2_fieldstrength;
%plots T2-decay curves for two tissues as well as their relative contrast
%http://en.wikipedia.org/wiki/Relaxation_(NMR)
%http://en.wikipedia.org/wiki/Spin-spin_relaxation_time

%Note T2 depends on tissue and field strength:
% Table 2.1 of Magnetic resonance imaging of the brain and spine, Volume 1 By Scott W. Atlas page 27
%T2 (ms)
%         Gray     White	CSF     Blood(Venous)   Blood(Arterial)   Fat
%   1.5T  100      80       2200    150             235                165 
%   3.0T  85       77       2200    80              165                133


x = 1:1:300;
xy15 = exp(-x./100);
xy30 = exp(-x./85);
p =plot(x,xy15,x,xy30);
xlabel('TE (ms)');
ylabel('Signal');
legend('1.5T','3.0T');
title( 'Gray matter T2 as a function of field strength');
set(p,'LineWidth',2);
set(gcf,'Color',[1 1 1]);
