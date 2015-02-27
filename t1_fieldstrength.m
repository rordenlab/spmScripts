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
%Alternatievely, see Wright et al. (2008) http://www.ncbi.nlm.nih.gov/pubmed/18259791
%   This work also lists measures from several prior studies
%         Gray    White Blood(a) Blood(v)
%   1.5T  1197     646  1412      1387
%   3.0T  1607     838  1613      1587
%   7.0T  1939    1126  2149      2120 
%Blood values from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3752414/
%  Blood (arterial) and venous T1(a) = 133.98 T + 1211.4, T1(v) = 133.27 T + 1187.4
%also nice
% http://www.vuiis.vanderbilt.edu/~welcheb/ExamCardDescriptions/tissue_nulling.html
%Alternatievely, see Rooney et al. (2007) http://www.ncbi.nlm.nih.gov/pubmed/17260370
%   This work also lists measures from several prior studies
%         cGray    White   CSF
%   1.5T  1188      656   4070
%   3.0T  1763      847
%   7.0T  2132     1200   4425 



x = 10:10:6000;
z15 = 1-exp(-x./1197);
z30 = 1-exp(-x./1607);
z70 = 1-exp(-x./1939);
p =plot(x,z15,x,z30,x,z70);
xlabel('TR (ms)');
ylabel('Recovery');
legend('Gray Matter (1.5T)','Gray Matter (3.0T)','Gray Matter (7.0T)');
title( 'Longitudinal Relaxation (T1) as function of field strength');
set(p,'LineWidth',2)
set(gcf,'Color',[1 1 1])



