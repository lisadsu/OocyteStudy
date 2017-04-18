% make paper figures
% Livia Zarnescu
% 1-27-15

%% Figure 1B

% load in aspiration depth data
load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse embryo\2-20-14 analysis\AutoMeasure\aspiration_data_2_20_14_E35.mat')
aspiration_depth = aspiration_depth * 40 * 10^-6 / 108;
xdata = t;
ydata = aspiration_depth;


% load in pressure data
fid = fopen('C:/Users/Livia/Desktop/IVF/Raw Data/Videos/Mouse Embryo/videos 2-20-14/E35_pressure.txt');
P = fscanf(fid,'%f');
P = P(43:end);
fclose(fid);
fid = fopen('C:/Users/Livia/Desktop/IVF/Raw Data/Videos/Mouse Embryo/videos 2-20-14/E35.txt');
tP = fscanf(fid, '%f');
tP = tP(43:end);
tP = cumsum(tP);
tP = tP / 1000;
tP = tP - min(tP);
fclose(fid);


% plot
x1 = xdata(2:end);
x2 = [-.1 0 xfine - min(xfine)];
y1 = 10^6*ydata(2:end);
y2 = [0 0 10^6*yfit];
x1 = [x1 NaN*ones(1, length(x2) - length(x1))];
y1 = [y1 NaN*ones(1, length(y2) - length(y1))];


[haxes, hline1, hline2] = plotyy([x1', x2'] , [y1', y2'], [-.1 0 tP'], [0 0 P']);
xlim(haxes(2), [-.05 .5]);
xlim(haxes(1), [-.05 .5]);
set(hline1(1), 'marker', 'o');
set(hline1(1), 'linestyle', 'none');
set(hline1(2), 'color', [0 0 1]);
set(haxes(2), 'ylim', [0 .8]);
set(haxes(1), 'ylim', [10 24]);
set(haxes(2), 'xtick', [])
set(haxes, 'fontsize', 14)
xlabel('time (seconds)')
ylabel(haxes(1), 'aspiration depth (\mum)')
ylabel(haxes(2), 'pressure applied (psi)');
title('Embryo Aspiration Depth into Micropipette')
set(haxes(2), 'ytick', [.1 .2 .3 .4 .5 .6 .7])
set(haxes(1), 'ytick', [12 14 16 18 20 22 24])
set(hline2, 'linewidth', 1)
set(hline1, 'linewidth', 1)
legend('Measured Depth', 'Fit to Model', 'Applied Pressure', 'Location', 'SouthEast')







%% 2B : Aspiration depth curves for mouse embryos


figure;
set(gca, 'FontSize', 14);
plot(xfine, 10^6*mean(gsum,1), 'Color', [0 .6 0], 'LineWidth', 2);
hold on;
plot(xfine, 10^6*mean(rsum,1), 'Color', [0 0 .6], 'LineWidth', 2);
h = legend('Blastocyst', 'No Blastocyst');
set(h, 'EdgeColor', [1 1 1]);
jbfill(xfine, 10^6*(mean(gsum,1)+std(gsum,[],1)), 10^6*(mean(gsum,1)-std(gsum,[],1)), ...
    [0 .6 0], 'none', [], .3);
jbfill(xfine, 10^6*(mean(rsum,1)+std(rsum,[],1)), 10^6*(mean(rsum,1)-std(rsum,[],1)), ...
    [0 0 .6], 'none', [], .3);

xlim([min(xfine) max(xfine)]);
xlabel('time (seconds)');
ylabel('aspiration depth (\mum)');
title('Average Aspiration Depth Curves');


% also plot a few specific examples
% 7-1-13 # 7
% 2-20-14 # 34

%% 2C : Scatterplot of mouse embryos



numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) & ...
    ~isnan(mList);% & ~isnan(timeLapseOut(:,1))';
mN = mList(numsToPlot);
k1N = k1list(numsToPlot);
% k1N = (k1N - min(k1N))/(max(k1N)-min(k1N));
n1N = n1list(numsToPlot);
% n1N = (n1N - min(n1N))/(max(n1N)-min(n1N));
tN = taulist(numsToPlot);
% tN = (tN - min(tN))/(max(tN)-min(tN));
k0N = k0list(numsToPlot);
% k0N = (k0N - min(k0N))/(max(k0N)-min(k0N));
tLN = timeLapseOut(numsToPlot,:);
n0N = n0list(numsToPlot);
eN = elonglist(numsToPlot);

T1 = tLN(:,1);
T2 = tLN(:,2);
T3 = tLN(:,3);

figure;
h = scatter3(k1N, n1N, k0N, 80, colorMat(numsToPlot,:), 'filled');
% h = scatter3(T1, T2, T3, 100, colorMat(numsToPlot,:), 'filled');

set(h, 'Marker', 'o');
set(gca, 'FontSize', 18);
title('3D scatter plot of parameters');
xlabel('k1 parameter');
ylabel('n1 parameter');
zlabel('k0 parameter');
% xlabel('cell cycle #1');
% ylabel('cell cycle #2');
% zlabel('cell cycle #3');
% axis([min(k1N) max(k1N) min(tN) max(tN)]);
axis([min(k1N) max(k1N) min(n1N) max(n1N) min(k0N) max(k0N)]);
% axis([min(T1) max(T1) min(T2) max(T2) min(T3) max(T3)]);
% axis([0 1 0 1 0 1]);
view(-23,8);
camlight right;
drawnow;
hold on;
set(gca, 'xscale', 'linear');
set(gca, 'yscale', 'log');
set(gca, 'zscale', 'log');

mVal = 5;
aAll = (1:length(k1list))';
a = aAll(numsToPlot & mList == mVal);% & mList == mVal);% & c2);
a = a - min(a) + 1;
b = num2str(a);
c = cellstr(b);
dx = -0.005; dy = 0.01; dz = .005; % displacement so the text does not overlay the data points


hold on;

aAll(mList == mVal)'

% ======== MAKE LEGEND =========
handleList = [];
labelList = cell(1,2);
colorList = [0 0 .6; ...    % dark blue
    0 .6 0];       % green
%
mV = [2 4];
for i = 1:2
    
    %     currHandle = scatter(k1N(find(mN == mV(i),1,'first')), ...
    %         n1N(find(mN == mV(i),1,'first')), ...
    %         100, colorList(i,:), 'filled');
    
    currHandle = scatter3(k1N(find(mN == mV(i),1,'first')), ...
        n1N(find(mN == mV(i),1,'first')), k0N(find(mN == mV(i),1,'first')), ...
        80, colorList(i,:), 'filled');
    
    %     currHandle = scatter3(T1(find(mN == mV(i),1,'first')), ...
    %         T2(find(mN == mV(i),1,'first')), T3(find(mN == mV(i),1,'first')),...
    %         100, colorList(i,:), 'filled');
    
    handleList = [handleList currHandle];
    
    if i == 1
        labelList{1} = 'no blastocyst';
    elseif i == 2
        labelList{2} = 'blastocyst';
    end
end

legend(handleList, labelList, 'Location', 'North');
legend('boxoff');
legend('boxon');
legend(handleList, labelList{1}, labelList{2}, ...
    'Location', 'North');


%% 2D: bar plot of live birth experiments

goodPercentage = [64 67 71 80];
controlPercentage = [43 50 50 53];
badPercentage = [21 17 17 40];

totalEmbryos = [14 12 14 15];
goodIndividual = [ones(1,9) zeros(1,5) ones(1,8) zeros(1,4) ones(1,10) zeros(1,4) ones(1,12) zeros(1,3)];
controlIndividual = [ones(1,6) zeros(1,8) ones(1,7) zeros(1,5) ones(1,6) zeros(1,8) ones(1,7) zeros(1,8)];
badIndividual = [ones(1,3) zeros(1,11) ones(1,2) zeros(1,10) ones(1,2) zeros(1,12) ones(1,6) zeros(1,9)];

x1 = [repmat('v', length(badIndividual),1); ...
     repmat('n', length(controlIndividual),1)];
x2 = [ones(sum(badIndividual),1); zeros(55 - sum(badIndividual), 1); ...
    ones(sum(controlIndividual),1); zeros(55 - sum(controlIndividual), 1)];
[tbl,chi2stat,pval] = crosstab(x1, x2)

figure;
h1 = bar(1, mean(badPercentage), 'facecolor', [.85 .65 .2]);
hold on;
e1 = terrorbar(1,mean(badPercentage),std(badPercentage),.1);
set(e1,'color', 'k', 'linewidth', 2);

h2 = bar(2, mean(controlPercentage), 'facecolor', [.6 .6 .6]);
e2 = terrorbar(2,mean(controlPercentage),std(controlPercentage),.1);
set(e2,'color', 'k', 'linewidth', 2);


h3 = bar(3, mean(goodPercentage), 'facecolor', [.2 .6 .9]);
e3 = terrorbar(3,mean(goodPercentage),std(goodPercentage),.1);
set(e3,'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [1 2 3])
set(gca, 'fontsize', 14);
ylim([0 100]);
set(gca, 'xticklabel', {'Non-viable', 'Control', 'Viable'});
xlabel('Predicted Embryo Viability');
ylabel('Embryos Resulting in Live Birth (%)');
title('Mechanical Parameters Predict Live Birth');
xlim([0.5 3.5]);
grid on;

%% 2E: Scatterplot of decBndDist vs percentage viable in live birth expts

goodDists = [2.3359 2.4707 0.8963 1.5148];
badDists = [4.0783 4.1445 4.4961 6.0862];

goodExpViable = [.8308 .8158 .9236 .9089];
badExpViable = [.7325 .6957 .6916 .7179];

goodPercentage = [64 67 71 80];
controlPercentage = [43 50 50 53];
badPercentage = [21 17 17 40];

figure, scatter([badExpViable goodExpViable], [badPercentage goodPercentage], ...
    'linewidth', 2);

p = polyfit([badExpViable goodExpViable], [badPercentage goodPercentage] , 1);
expViableFit = polyval(p, [0 badExpViable goodExpViable 1]);

[xValSort, I] = sort([0 badExpViable goodExpViable 1], 'ascend');

hold on;
plot(xValSort, expViableFit(I), 'linewidth', 2);

% rsq = .8676
rsq = 1 - sum((expViableFit - [0 badPercentage goodPercentage 1]).^2) / ...
    sum(([0 badPercentage goodPercentage 1] - mean([0 badPercentage goodPercentage 1])).^2);

axis([.6 1 0 100]);
set(gca, 'fontsize', 14);
xlabel('expected % forming blastocysts');
ylabel('actual % resulting in live birth');
grid on;
title('Viability Predictions Correlate with Live Birth Rates');



%% Figure 3A

fileName = 'C:/Users/Livia/Desktop/IVF/Raw Data/Videos/Human/human tests 3-19-14/E';

for i = 1:20
    
    A = imread([fileName num2str(i) '.jpg']);
    
    if length(size(A)) == 3
        imwrite(medfilt2(rgb2gray(A), [3 3]), [fileName num2str(i) '.jpg']);
    end
    
end

%% Figure 4C/D (ROC and PR curves for human, mouse, mech, cell cycle feature selection)

clear all;
load('forwardFeatureSelection.mat');
Xref = 0:.01:1;

close all
figure(1); 
set(gca, 'fontsize', 14);
% h1 = plot(ROCy_mouse_mech, PRy_mouse_mech, 'linewidth', 2);
hold on;
h2 = plot(ROCy_human_mech, PRy_human_mech, 'linewidth', 2, 'linestyle', '--');
% h3 = plot(ROCy_mouse_all, PRy_mouse_all, 'linewidth', 2, 'color', [1 0 0]);
h4 = plot(ROCy_human_all, PRy_human_all, 'linewidth', 2, 'color', [1 0 0], 'linestyle', '--');
% h5 = plot(ROCy_mouse_cellCycle, PRy_mouse_cellCycle, 'linewidth', 2, 'color', [0 .6 0]);
h6 = plot(ROCy_human_cellCycle, PRy_human_cellCycle, 'linewidth', 2, 'color', [0 .6 0], 'linestyle', '--');
% legend([h1 h5 h3 h2 h6 h4], 'Mouse, mech only', 'Mouse, cell cycle only', 'Mouse, all params', ...
%      'Human, mech only', 'Human, cell cycle only','Human, all params', 'Location', 'SouthWest');
 legend([h2 h6 h4], 'Human, mech only', 'Human, cell cycle only', ...
     'Human, all params', 'Location', 'SouthWest');
axis([0 1 0 1]);
grid on;
xlabel('Recall (sensitivity)');
ylabel('Precision (PPV)');
title('Precision-Recall curves');

figure(2); 
set(gca, 'fontsize', 14);
% h7 = plot(Xref, ROCy_mouse_mech, 'linewidth', 2);
hold on;
h8 = plot(Xref, ROCy_human_mech, 'linewidth', 2, 'linestyle', '--');
% h9 = plot(Xref, ROCy_mouse_all, 'linewidth', 2, 'color', [1 0 0]);
h10 = plot(Xref, ROCy_human_all, 'linewidth', 2, 'color', [1 0 0], 'linestyle', '--');
% h11 = plot(Xref, ROCy_mouse_cellCycle, 'linewidth', 2, 'color', [0 .6 0]);
h12 = plot(Xref, ROCy_human_cellCycle, 'linewidth', 2, 'color', [0 .6 0], 'linestyle', '--');
% legend([h7 h11 h9 h8 h12 h10], 'Mouse, mech only', 'Mouse, cell cycle only', 'Mouse, all params', ...
%      'Human, mech only', 'Human, cell cycle only','Human, all params', 'Location', 'SouthEast');
legend([h8 h12 h10], 'Human, mech only', 'Human, cell cycle only', ...
    'Human, all params', 'Location', 'SouthEast');
axis([0 1 0 1]);
grid on;
xlabel('1 - Specificity');
ylabel('Sensitivity');
title('ROC curves');

%% Figure 4E/F: Forward feature selection

clear all;
close all;
load('forwardFeatureSelection.mat');
Xref = 0:.01:1;

figure(1); 
set(gca, 'fontsize', 14);
% h1 = plot(ROClist_mouse_mech, 'linewidth', 2);
hold on;
% h2 = plot(PRlist_mouse_mech, 'linewidth', 2, 'color', [1 0 0]);
h3 = plot(ROClist_human_mech, 'linewidth', 2, 'linestyle', '--');
h4 = plot(PRlist_human_mech, 'linewidth', 2, 'linestyle', '--', 'color', [1 0 0]);
% legend([h1 h2 h3 h4], 'Mouse, AUC_R_O_C', 'Mouse, AUC_P_R', ...
%     'Human, AUC_R_O_C', 'Human, AUC_P_R', 'Location', 'SouthWest');
legend([h3 h4], 'Human, AUC_R_O_C', 'Human, AUC_P_R', 'Location', 'SouthWest');
xlabel('feature number');
ylabel('Area Under Curve');
title('Feature Selection, Mechanical Parameters');
ylim([.3 1]);
set(gca, 'xtick', [1 2 3 4]);
grid on;
box on;

figure(2); 
set(gca, 'fontsize', 14);
h1 = plot(ROClist_mouse_all, 'linewidth', 2);
hold on;
h2 = plot(PRlist_mouse_all, 'linewidth', 2, 'color', [1 0 0]);
h3 = plot(ROClist_human_all, 'linewidth', 2, 'linestyle', '--');
h4 = plot(PRlist_human_all, 'linewidth', 2, 'linestyle', '--', 'color', [1 0 0]);
legend([h1 h2 h3 h4], 'Mouse, AUC_R_O_C', 'Mouse, AUC_P_R', ...
    'Human, AUC_R_O_C', 'Human, AUC_P_R', 'Location', 'SouthWest');
xlabel('feature number');
ylabel('area under curve');
title('Feature Selection, All Parameters');
ylim([.6 1]);
xlim([1 7]);
set(gca, 'xtick', [1 2 3 4 5 6 7]);
grid on;


%% Figure 4G: AUC_ROC and AUC_PR over time

clear all;
% close all;
load('forwardFeatureSelection.mat');

numCells = [1 2 3 4];
% AUC_ROC_mouse = [max(ROClist_mouse_mech) max(ROClist_mouse_5feat) ...
%     max(ROClist_mouse_6feat) max(ROClist_mouse_all)];
% AUC_PR_mouse = [max(PRlist_mouse_mech) max(PRlist_mouse_5feat) ...
%     max(PRlist_mouse_6feat) max(PRlist_mouse_all)];

AUC_ROC_human_all = [max(ROClist_human_mech) max(ROClist_human_5feat) ...
    max(ROClist_human_6feat) max(ROClist_human_all)];
% AUC_PR_human = [max(PRlist_human_mech) max(PRlist_human_5feat) ...
%     max(PRlist_human_6feat) max(PRlist_human_all)];

AUC_ROC_human_mech = max(ROClist_human_mech)*ones(1,4);
AUC_ROC_human_cc = [0 0.72 0.92 0.95];


figure(1); 
clf;
set(gca, 'fontsize', 14);
% h1 = plot(numCells, AUC_ROC_mouse, 'linewidth', 2);
hold on;
% h2 = plot(numCells, AUC_PR_mouse, 'linewidth', 2, 'color', [1 0 0]);
h3 = plot(numCells, AUC_ROC_human_all, 'linewidth', 2, 'linestyle', '--', 'color', [1 0 0]);
h4 = plot(numCells, AUC_ROC_human_mech, 'linewidth', 2, 'linestyle', '--');
h5 = plot(numCells, AUC_ROC_human_cc, 'linewidth', 2, 'linestyle', '--', 'color', [0 .6 0]);
% h4 = plot(numCells, AUC_PR_human, 'linewidth', 2, 'linestyle', '--', 'color', [1 0 0]);
% legend([h1 h2 h3 h4], 'Mouse, AUC_R_O_C', 'Mouse, AUC_P_R', ...
%     'Human, AUC_R_O_C', 'Human, AUC_P_R', 'Location', 'SouthEast');
legend([h4 h5 h3], 'mechanics only', 'cell cycle only', ...
    'mechanics + cell cycle', 'Location', 'SouthEast');
xlabel('number of cells in embryo');
ylabel('area under ROC curve');
title('Predictive Power Increases Over Time');
ylim([0 1]);
set(gca, 'ytick', [0 .2 .5 .7 .8 .9 1.0]);
set(gca, 'yticklabel', {'0', '', '', '0.7', '0.8', '0.9', '1'});
set(gca, 'xtick', [1 2 3 4]);
grid on;
box on;

%% Figure Bar plot of parameters +/- stdev

colorList = [0 0 .6; ...     % dark blue
             0 .6 .0; ...
             0 0 .6; ...
             0 .6 0];       % green

         
figure;
k11 = bar(.75, mean(k1N(mN == 1)), .4, 'facecolor', colorList(1,:));
hold on;
ek11 = errorbar(.75, mean(k1N(mN == 1)),std(k1N(mN == 1)),'color', 'k', 'linewidth', 2);
k12 = bar(1.25, mean(k1N(mN == 4)), .4, 'facecolor', colorList(2,:));
ek12 = errorbar(1.25, mean(k1N(mN == 4)),std(k1N(mN == 4)),'color', 'k', 'linewidth', 2);

n11 = bar(1.75, mean(log(n1N(mN == 1))), .4, 'facecolor', colorList(1,:));
en11 = errorbar(1.75, mean(log(n1N(mN == 1))),std(log(n1N(mN == 1))),'color', 'k', 'linewidth', 2);
n12 = bar(2.25, mean(log(n1N(mN == 4))), .4, 'facecolor', colorList(2,:));
en12 = errorbar(2.25, mean(log(n1N(mN == 4))),std(log(n1N(mN == 4))),'color', 'k', 'linewidth', 2);

k01 = bar(2.75, mean(log(k0N(mN == 1))), .4, 'facecolor', colorList(1,:));
ek01 = errorbar(2.75, mean(log(k0N(mN == 1))),std(log(k0N(mN == 1))),'color', 'k', 'linewidth', 2);
k02 = bar(3.25, mean(log(k0N(mN == 4))), .4, 'facecolor', colorList(2,:));
ek02 = errorbar(3.25, mean(log(k0N(mN == 4))),std(log(k0N(mN == 4))),'color', 'k', 'linewidth', 2);

t1 = bar(3.75, mean(tN(mN == 1)), .4, 'facecolor', colorList(1,:));
et1 = errorbar(3.75, mean(tN(mN == 1)),std(tN(mN == 1)),'color', 'k', 'linewidth', 2);
t2 = bar(4.25, mean(tN(mN == 4)), .4, 'facecolor', colorList(2,:));
et2 = errorbar(4.25, mean(tN(mN == 4)),std(tN(mN == 4)),'color', 'k', 'linewidth', 2);

% set(gca, 'xtick', [1 2 3])
% set(gca, 'fontsize', 14);
% ylim([0 100]);
% set(gca, 'xticklabel', {'Non-viable', 'Control', 'Viable'});
% xlabel('Predicted Embryo Viability');
% ylabel('Embryos Resulting in Live Birth (%)');
% title('Mechanical Parameters Predict Live Birth');
% xlim([0.5 3.5]);
%          


%% Spectrum presentation

% make smoothed kernel density plot
[fN xiN] = ksdensity(decDist(mN == 1));
[fV xiV] = ksdensity(decDist(mN == 4));

figure;
set(gca, 'FontSize', 14);
title('Embryo Viability Likelihood');
xlabel('decision boundary distance');
ylabel('histogram density');

hold on;
v1 = plot(xiV, fV, 'Color', [.3 .8 .3], 'linewidth', 2);
nv1 = plot(xiN, fN, 'Color', [.8 .3 .3], 'linewidth', 2);
h = area(xiV, fV, 'EdgeColor', [.3 .8 .3], 'FaceColor', [.3 .8 .3]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .3);
h = area(xiN, fN, 'EdgeColor', [.8 .3 .3], 'FaceColor', [.8 .3 .3]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .3);
grid on;
box on;

%% Control viability

controlV = zeros(1,35);
controlV(1:23) = 100;
controlSE = 196*sqrt((23/35)*(12/35)/35); % 95% CI

testV = zeros(1,282);
testV(1:197) = 100;
testSE = 196*sqrt((197/282)*(85/282)/282); % 95% CI


[h p] = ttest2(controlV, testV)

figure;
h1 = bar(1, mean(controlV), 'facecolor', [.5 .5 .5]);
hold on;
e1 = terrorbar(1, mean(controlV), 2*controlSE, .1);
set(e1, 'color', 'k', 'linewidth', 2);


h2 = bar(2, mean(testV), 'facecolor', [.5 .5 .5]);
e2 = terrorbar(2,mean(testV),2*testSE, .1);
set(e2, 'color', 'k', 'linewidth', 2);



set(gca, 'xtick', [1 2])
set(gca, 'fontsize', 14);
set(gca, 'xticklabel', {'Control', 'Measured'});
ylabel('% blastocyst formation');
title('Mechanical measurement does not affect viability');
xlim([0.5 2.5]);
ylim([0 100]);
grid on;

%% generate sample histogram for Figure 1

a1 = normrnd(1.5,0.3,1,100);
a2 = normrnd(0.8,0.2,1,100);

% make smoothed kernel density plot
[fN xiN] = ksdensity(a1);
[fV xiV] = ksdensity(a2);

figure;
hold on;
v1 = plot(xiV, fV, 'Color', [.8 .6 .3], 'linewidth', 2);
nv1 = plot(xiN, fN, 'Color', [.3 .6 .8], 'linewidth', 2);
h = area(xiV, fV, 'EdgeColor', [.8 .6 .3], 'FaceColor', [.8 .6 .3]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .3);
h = area(xiN, fN, 'EdgeColor', [.3 .6 .8], 'FaceColor', [.3 .6 .8]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .3);
xlim([-.1,2.5]);
grid on;
box on;


%% Other figures...


% Fig 4C: makeMicroinjectionFigs.m
% Fig 4D: analyzeMicroinjectionData.m
% Fig 4E: makeMicroinjectionFigs.m
% Fig 4F: plotAllMouseOocytesSoFar.m



%% Fig 3E and F

% reregenerate high-res IPA fig outputs

downColor = [0 1 0];
upColor = [1 0 0];

panelE = {'BRCA1', 'RAE1', 'MBD4', 'FOXO3', 'BUB3', 'MRE11A', 'MCM7', ...
    'BUB1', 'TTI1', 'FANCC', 'RAD18', 'NBN', 'SMC3', 'TERF1', 'HUS1', ...
    'RBBP8', 'DBF4', 'NDC80', 'AURKB', 'NABP2', 'ATR', 'UIMC1', 'PARP2', ...
    'ZW10', 'CDK2', 'XPC', 'CHEK2', 'TP53BP1', 'PIM1', 'BUB1B', 'DDB1', ...
    'ZWILCH', 'FANCD2', 'CDKN1A', 'CDC6', 'HELB', 'E2F1', 'PTTG1', 'ZWINT', ...
    'DCLRE1B', 'ZAK', 'SIN3B', 'HSPA1A', 'HSPA1B', 'CKS2', 'TP53', 'CCNG2', ...
    'MXD1', 'MAPK14', 'PLK3', 'MAD2L1', 'CCNA2', 'NCAPD3', 'NCAPD2', ...
    'CAPN2', 'SMC4', 'DLGAP5', 'NCAPG', 'HCFC1', 'KIF2C', 'KIF15', ...
    'KIF18A', 'KIFC1'};
logFCE = zeros(1,length(panelE));
colorE = zeros(length(panelE),3);

panelF = {'SURF4', 'PWWP2A', 'NEK6', 'DEC1', 'BEND3', 'MIER3', 'MIER2', 'HDAC1', ...
    'ZBTB2', 'ZNF609', 'TRERF1', 'ARID5B', 'MT1E', 'SRD5A1', 'ZNF462', 'ZNF620', ...
    'CBX5', 'ZNF45', 'MECP2', 'ENPP4', 'DLX6', 'MKI67', 'TSPAN14', 'REST', ...
    'GON4L', 'YY1', 'ACTR8', 'INO80D', 'ZNF148', 'HDAC3', 'PARP14', 'PRDM4', ...
    'WDTC1', 'MXD1' 'RNF17'};
logFCF = zeros(1, length(panelF));
colorF = zeros(length(panelF), 3);

fileE = 'C:\Users\Livia\Dropbox\Embryo Mechanics outline shared\Data\IPA\logFC.txt';
f = fopen(fileE, 'r');
c = textscan(f, '%s %f');
fclose(f);

for i = 1:length(c{1})
    nameCurr = c{1}(i);
    for j = 1:length(panelE)
        if isequal(nameCurr{1}(2:end-1), panelE{j})
            logFCE(j) = c{2}(i);
        end
    end
end


for i = 1:length(logFCE)
    if logFCE(i) > 0
        colorE(i,:) = [1 1 1] - [1 0 1]*logFCE(i);
    else
        colorE(i,:) = [1 1 1] + [0 1 1]*logFCE(i);
    end
    colorE(colorE < 0) = 0;
end


for i = 1:length(c{1})
    nameCurr = c{1}(i);
    for j = 1:length(panelF)
        if isequal(nameCurr{1}(2:end-1), panelF{j})
            logFCF(j) = c{2}(i);
        end
    end
end


for i = 1:length(logFCF)
    if logFCF(i) > 0
        colorF(i,:) = [1 1 1] - [1 0 1]*logFCF(i);
    else
        colorF(i,:) = [1 1 1] + [0 1 1]*logFCF(i);
    end
    colorF(colorF < 0) = 0;
end


%% now plot
figure(1);
clf;

for i = 1:length(logFCE)
   rectangle('position', [80*mod(i,8) 4*floor(i/8)  8*length(panelE{i}) 2], ...
       'edgecolor', 'k', 'linewidth', 1, 'facecolor', colorE(i,:), ...
       'curvature', 1);
   hold on;
   text(80*mod(i,8)+1.5*length(panelE{i}), 4*floor(i/8)+1, panelE{i}, 'fontsize', 8);
end
    
axis([-50 650 -5 35])


figure(2);
clf;

for i = 1:length(logFCF)
   rectangle('position', [80*mod(i,8) 4*floor(i/8)  8*length(panelF{i}) 2], ...
       'edgecolor', 'k', 'linewidth', 1, 'facecolor', colorF(i,:), ...
       'curvature', 1);
   hold on;
   text(80*mod(i,8)+1.5*length(panelF{i}), 4*floor(i/8)+1, panelF{i}, 'fontsize', 8);
end
    
axis([-50 650 -5 35])











