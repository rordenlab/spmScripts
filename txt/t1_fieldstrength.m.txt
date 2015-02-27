function t1_fieldstrength;
%plots T1-recovery for same tissue at two field strengths
%http://en.wikipedia.org/wiki/Spin-lattice_relaxation_time
%
%Note T1 depends on tissue and field strength:
%  Magnetic resonance imaging of the brain and spine, Volume 1 By Scott W. Atlas page 27
%  T1 (ms)
%         Gray    White   CSF  Blood  Fat
%   1.5T  1000     710   4000  1435   300
%   3.0T  1331     832   4000  1584   380


x = 10:10:6000;
z15 = 1-exp(-x./1000);
z30 = 1-exp(-x./1331);
p =plot(x,z15,x,z30);
xlabel('TR (ms)');
ylabel('Recovery');
legend('Gray Matter (1.5T)','Gray Matter (3.0T)');
title( 'Longitudinal Relaxation (T1) as function of field strength');
set(p,'LineWidth',2)
set(gcf,'Color',[1 1 1])
