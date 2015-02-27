function t2_tissue;
%plots T2-decay curves for several tissues at the same field strength
%http://en.wikipedia.org/wiki/Relaxation_(NMR)
%http://en.wikipedia.org/wiki/Spin-spin_relaxation_time

%Note T2 depends on tissue and field strength:
% Table 2.1 of Magnetic resonance imaging of the brain and spine, Volume 1 By Scott W. Atlas page 27

%T2 (ms)
%         Gray     White	CSF     Blood(Venous)   Blood(Arterial)   Fat
%   1.5T  100      80       2200    150             235                165 
%   3.0T  85       77       2200    80              165                133


x = 1:1:500;
xyGM = exp(-x./85);
xyWM = exp(-x./77);
xyBlood = exp(-x./130);
xyCSF = exp(-x./2200);
p =plot(x,xyGM,x,xyWM,x,xyBlood,x,xyCSF);
xlabel('TE (ms)');
ylabel('Signal');
legend('Gray Matter','White Matter','Blood','CSF');
title( 'Dephasing: Transverse Relaxation (T2) at 3.0T');
set(p,'LineWidth',2);
set(gcf,'Color',[1 1 1]);
