% Make figures for review paper
% Livia Zarnescu Yanez
% 3-30-16
%
% Compare values of oocyte and embryo Young's moduli in literature


% 1 = mouse, 2 = human, 3 = bovine, 4 = hamster
xvals = [1.0 1.1 1.2 1.3 1.4 1.5 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7];
species = [1 1 1 2 3 4 1 1 1 1 1 2 3 4]';
% colorVals = [35 50 159; 231 174 24; 16 155 106; 231 91 24]/255;
colorVals = [50 50 150; 150 150 200; 180 50 40; 220 140 140]/255;
    
% softening ratio, immature to mature oocytes
oocyteMaturationStiffnessRatio = [3.6/4.2, 22.8/8.26, .0844/.0563, (7.2+9.3)/(3.1+11.2), 89/22, 6.2/2.5];
oocyteMaturationStiffnessRatioRef = ['Drobnis1988', 'Murayama2006', 'Yanez2016', 'Abadie2014', 'Papi2010', 'Drobnis1988'];

% % softening ratio, high to low quality oocytes
% oocyteQualityStiffnessRatio = [3.1/1.6];
% oocyteQualityStiffnessRatioRef = ['Liu2010'];

% hardening ratio, oocyte to embryo
zonaHardeningRatio = [4.1/2.3, 42.2/17.9, 22.3/8.26, 36.9/11.8, 2.8, (13.18+14.19)/(7.34+7.47), 84/22, 5.5/4.95];
zonaHardeningRatioRef = ['Drobnis1988', 'Sun2003', 'Murayama2006', 'Khalilian2010', 'Yanez2016', 'Khalilian2010', 'Papi2010', 'Drobnis1988'];

allVals = [oocyteMaturationStiffnessRatio zonaHardeningRatio];

% 1. Plot maturation softening


figure(1); clf;
hold on;
set(gca, 'fontsize', 14);

for i = 1:length(species)
    h(i) = bar(xvals(i), allVals(i), .08, 'facecolor', colorVals(species(i),:));
end

xlim([.8 2.9]);
ylim([0 6]);
grid on;
set(gca, 'xticklabel', '');
% set(gca, 'xtick', [1.25 2.3]);
% set(gca, 'xticklabel', {'oocyte maturation', 'zona hardening'})
