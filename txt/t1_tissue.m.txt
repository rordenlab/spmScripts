function t1_tissue;
%plots T1-recovery for two tissues at the same field strength
%http://en.wikipedia.org/wiki/Spin-lattice_relaxation_time
%
%Note T1 depends on tissue and field strength:
%  Magnetic resonance imaging of the brain and spine, Volume 1 By Scott W. Atlas page 27
%  T1 (ms)
%         Gray    White   CSF  Blood  Fat
%   1.5T  1000     710   4000  1435   300
%   3.0T  1331     832   4000  1584   380

x = 10:10:6000;
zGM = 1-exp(-x./1331);
zWM = 1-exp(-x./832);
p =plot(x,zGM,x,zWM);
xlabel('TR (ms)');
ylabel('Recovery');
legend('Gray Matter','White Matter');
title( 'Longitudinal Relaxation (T1) at 3.0T');
set(p,'LineWidth',2)
set(gcf,'Color',[1 1 1])
