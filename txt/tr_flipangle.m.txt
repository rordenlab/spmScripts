function tr_flipangle(tr);
%Plots optimal flip angle as a function of repetition time (TR) for given T1
%
%http://en.wikipedia.org/wiki/Ernst_angle
%http://desoft03.usc.es/rmnweb/calibrations/Ernstcalculator.htm
%From  Mark Cohen's http://airto.ccn.ucla.edu/BMCweb/HowTo/CalcFlip.html
% ported to Matlab by Chris Rorden
%
%Very short TRs (relative to the tissue's T1) will lead to decreased signal 
% because the recovery T1 recovery will be incomplete. In these case, a
% shallower flip angle will offer MORE signal, as it prevents saturation.
% This code estimates the optimal flip angle across TRs for a given T1.
%
%Note T1 depends on tissue and field strenght:
%  Magnetic resonance imaging of the brain and spine, Volume 1 By Scott W. Atlas page 27
%  T1 (ms)
%         Gray    White   CSF  Blood  Fat
%   1.5T  1000     710   4000  1435   300
%   3.0T  1331     832   4000  1584   380
%For fMRI studies please read http://www.ncbi.nlm.nih.gov/pubmed/21073963
T1 = 1400; %assumed T1 for brain tissue


if nargin <1 %no Positive files
    tr = 1550; %you can set the tr to any value you wish
end;

%compute a single optimal flip angle
fprintf('Assuming a brain T1 of %.2fms and a TR of %.2fms, the optimum flip angle is %.2f and the signal is %.0f%% relative to infinite TR.\n',T1, tr,ErnstFlipDeg(T1,tr),ErnstGain(T1,tr));

x = 10:10:3000;
y = ErnstFlipDeg(T1,x);
y2 = ErnstGain(T1,x);
[Ax,H1,H2]  = plotyy(x,y,x,y2);
set(get(Ax(1),'Ylabel'),'String','Optimal flip angle (deg)') 
set(get(Ax(2),'Ylabel'),'String','Signal relative to infinite TR (%)') 
xlabel('TR (ms)');
%legend('Flip angle (deg)','Signal relative to infinite TR (%)');
title( ['Ernst angle for T1 of ',int2str(T1), 'ms']);
set(H1,'LineWidth',2);
set(H2,'LineWidth',2);
set(gcf,'Color',[1 1 1])

function [flipAngle] = ErnstFlipDeg(T1,tr);
flipAngle =  acos( exp(-tr / T1));
flipAngle = 180 * flipAngle / pi; %radians to degrees

function [rSNR] = ErnstGain(T1,tr);
flipAngle =  acos( exp(-tr / T1));
%fprintf('Assuming a brain T1 of %.2fms and a TR of %.2fms, the optimum flip angle is %.2f and the gain is %.0f\n',T1, tr,ErnstFlipDeg(T1,tr),ErnstGain(T1,tr));
Ratio = sin( flipAngle ) / ( 1-cos(flipAngle) .* exp( -tr/T1 ) );
%fprintf('SNR improvement by using optimal flip angle: %.2f\n', 100.0 * (Ratio - 1) );
rSNR = 100 * Ratio * (1-exp( -tr/T1 ));
%fprintf('Relative signal (compared to infinite tr): %.2f, therefore SNR gain is %.0f\n', rSNR, 100-rSNR  );
%gain = 100-rSNR;

