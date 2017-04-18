% Make figures from microinjection experiments from 9-24-14
% incorporates both mechanics and cortical granule data
% Livia Zarnescu
% 1-27-15

%% START HERE FOR FIGS Make scatterplot from params

fileToSave = 'C:/Users/Livia/Desktop/IVF/Processed Data/Mouse Embryo/microinjection_data/embryoInfoMicroinjection.mat';
load(fileToSave);
circMeanList = embryoInfo.circMeanList;
circStdList = embryoInfo.circStdList;
circDiffList = embryoInfo.circDiffList;
radMeanList = embryoInfo.radMeanList;
brightParamList = embryoInfo.brightParamList;
ROImeanList = embryoInfo.ROImeanList;
groupList = embryoInfo.groupList;
embryosInGroup = embryoInfo.embryosInGroup;

groupsToPlot = [1 1 1 1 1 1 1 1 1 1 1 1];
groupNums = 1:12;

colorVec = [1 .6 .6; ... % Ab-soft: pale red
    .8 .3 .3; ... % Ab-mid: med red
    .6 0 0; ... % Ab-stiff: dark red
    .6 .6 1; ... % BAPTA-soft: pale blue
    .3 .3 .8; ... % BAPTA-mid: mid blue
    0 0 .6; ... % BAPTA-stiff: dark blue
    .6 1 .6; ... % Control-soft: pale green
    .3 .8 .3; ... % Control-mid: mid green
    0 .6 0; ... % Control-stiff: dark green
    1 1 .6; ... % NonControl-soft: pale yellow
    .8 .8 .3; ... % NonControl-med: mid yellow
    .6 .6 0]; % NonControl-stiff: dark yellow

colorsToPlot = [repmat(colorVec(1,:), embryosInGroup(1)*groupsToPlot(1), 1); ...
                repmat(colorVec(2,:), embryosInGroup(2)*groupsToPlot(2), 1); ...
                repmat(colorVec(3,:), embryosInGroup(3)*groupsToPlot(3), 1); ...
                repmat(colorVec(4,:), embryosInGroup(4)*groupsToPlot(4), 1); ...
                repmat(colorVec(5,:), embryosInGroup(5)*groupsToPlot(5), 1); ...
                repmat(colorVec(6,:), embryosInGroup(6)*groupsToPlot(6), 1); ...
                repmat(colorVec(7,:), embryosInGroup(7)*groupsToPlot(7), 1); ...
                repmat(colorVec(8,:), embryosInGroup(8)*groupsToPlot(8), 1); ...
                repmat(colorVec(9,:), embryosInGroup(9)*groupsToPlot(9), 1); ...
                repmat(colorVec(10,:), embryosInGroup(10)*groupsToPlot(10), 1); ...
                repmat(colorVec(11,:), embryosInGroup(11)*groupsToPlot(11), 1); ...
                repmat(colorVec(12,:), embryosInGroup(12)*groupsToPlot(12), 1)];
                    
figure, scatter3([circMeanList{groupNums(groupsToPlot == 1)}], ...
    [radMeanList{groupNums(groupsToPlot == 1)}], ...
    [brightParamList{groupNums(groupsToPlot == 1)}], 100, colorsToPlot, 'filled');
xlabel('x');
ylabel('y');
zlabel('z');


%% Load mech params and compute average 

date1 = '9-24-14';
date2 = '9_24_14';

embryoDataDirectory = ['C:\Users\Livia\Desktop\IVF\Processed Data\' ...
    'Mouse embryo\', date1, ' analysis\'];
dataFileName = ['aspiration_data_', date2, '_E'];
fullDataDir = [];

k0list = cell(1,12);
k1list = cell(1,12);
n1list = cell(1,12);

embryosInGroup = [4 7 5 2 8 5 5 5 6 2 6 7];
embryoNumsInGroups = {[15 28 1 14], [3 8 10 12 32 20 16], [11 18 34 27 19], ...
    [38 39], [36 40 41 43 53 54 56 57], [49 37 45 50 44], ...
    [68 71 75 72 84], [67 70 81 82 83], [78 69 80 66 79 77], ...
    [87 95], [85 86 93 96 100 104], [91 88 92 94 98 99 101]};

for i = 1:length(embryosInGroup)
    
    currEmbryoNums = [embryoNumsInGroups{i}];
    
    currk0 = [];
    currk1 = [];
    currn1 = [];
    
    for j = 1:length(currEmbryoNums)
        
        dataFileCurr = [embryoDataDirectory 'AutoMeasure\' dataFileName ...
                num2str(currEmbryoNums(j)) '.mat'];
        
        if ~exist(dataFileCurr, 'file')
            dataFileCurr = [embryoDataDirectory dataFileName ...
                num2str(currEmbryoNums(j)) '.mat'];
        end
        
        load(dataFileCurr);
        currk0 = [currk0 k0];
        currk1 = [currk1 k1];
        currn1 = [currn1 n1];
        
    end
    
    k0list{i} = currk0;
    k1list{i} = currk1;
    n1list{i} = currn1;
    
end


%% Make scatterplot between average signal intensity and average stiffness

groupsToPlot = [1 0 1 0 0 0 1 1 1 0 0 0];
groupNums = 1:12;

colorsToPlot = [repmat(colorVec(1,:), embryosInGroup(1)*groupsToPlot(1), 1); ...
                repmat(colorVec(2,:), embryosInGroup(2)*groupsToPlot(2), 1); ...
                repmat(colorVec(3,:), embryosInGroup(3)*groupsToPlot(3), 1); ...
                repmat(colorVec(4,:), embryosInGroup(4)*groupsToPlot(4), 1); ...
                repmat(colorVec(5,:), embryosInGroup(5)*groupsToPlot(5), 1); ...
                repmat(colorVec(6,:), embryosInGroup(6)*groupsToPlot(6), 1); ...
                repmat(colorVec(7,:), embryosInGroup(7)*groupsToPlot(7), 1); ...
                repmat(colorVec(8,:), embryosInGroup(8)*groupsToPlot(8), 1); ...
                repmat(colorVec(9,:), embryosInGroup(9)*groupsToPlot(9), 1); ...
                repmat(colorVec(10,:), embryosInGroup(10)*groupsToPlot(10), 1); ...
                repmat(colorVec(11,:), embryosInGroup(11)*groupsToPlot(11), 1); ...
                repmat(colorVec(12,:), embryosInGroup(12)*groupsToPlot(12), 1)];

circToPlot = zeros(1,12);
stiffToPlot = zeros(1,12);
figure;
hold on;
h = [];
plotSoFar = 0;

newColorVec = [.8 .6 .3; ...
               .8 .6 .3; ...
               .3 .3 .8; ...
               .3 .3 .8; ...
               .3 .3 .8];

for i = 1:12
    
    if (groupsToPlot(i))
        
        plotSoFar = plotSoFar + 1;
        
        circToPlot(i) = mean([ROImeanList{i}]);
        stiffToPlot(i) = mean([k1list{i}]);
        
        if i == 1
            stiffToPlot(i) = stiffToPlot(i) + .001; % for visibility
        end
        
        he = herrorbar(stiffToPlot(i), circToPlot(i), std([k1list{i}]));
        set(he, 'color', [0 0 0], 'linewidth', 2);
        
        he = line([stiffToPlot(i) stiffToPlot(i)], ...
            [circToPlot(i) + std([brightParamList{i}]) ...
            circToPlot(i) - std([brightParamList{i}])], ...
            'color', [0 0 0], 'linewidth', 2);
        
        % soft group: filled square
        % stiff group: filled circle
        % med: not filled
        if i == 1
            h(i) = scatter(stiffToPlot(i), circToPlot(i), 350, [.8 .6 .3], 'filled');
            set(h(i), 'Marker', 'square');
        elseif i == 3
            h(i) = scatter(stiffToPlot(i), circToPlot(i), 350, [.8 .6 .3]);
            set(h(i), 'Marker', 'diamond', 'LineWidth', 2);
        elseif i == 7
            h(i) = scatter(stiffToPlot(i), circToPlot(i), 350, [.3 .3 .8], 'filled');
            set(h(i), 'Marker', 'square');

        elseif i == 9
            h(i) = scatter(stiffToPlot(i), circToPlot(i), 350, [.3 .3 .8]);
            set(h(i), 'Marker', 'diamond', 'LineWidth', 2);
        
        elseif i == 8
            h(i) = scatter(stiffToPlot(i), circToPlot(i), 350, [.7 .7 .7], 'filled');
            set(h(i), 'LineWidth', 2);        
%         h(i) = scatter(stiffToPlot(i), circToPlot(i), 150, newColorVec(plotSoFar,:), 'filled');
        end
    end
end
      

set(gca, 'fontsize', 14);
ylabel('CG brightness');
xlabel('stiffness');
legend(h([1 3 7 9 8]), 'Antibody-soft', 'Antibody-stiff', 'Control-soft', ...
    'Control-stiff', 'Control-mid', 'Location', 'NorthWest');
grid on;
xlim([.1 .23]);

% stats tests
% 1. CG brightness btwn stiffest quartile of Ab-injected vs control-injected
[p h] = ranksum(circMeanList{groupNums == 3}, circMeanList{groupNums == 9})
% 2. Stiffness btwn stiffest quartile of Ab-injected vs control-injected
[p h] = ranksum(k1list{groupNums == 3}, k1list{groupNums == 9})
% 2. Stiffness btwn softest quartile of Ab-injected vs control-injected
[p h] = ranksum(k1list{groupNums == 1}, k1list{groupNums == 7})
% 3. CG brightness btwn stiffest quartile of Ab-injected vs middle quartile of control-injected
[p h] = ranksum(circMeanList{groupNums == 3}, circMeanList{groupNums == 8})
% 3. CG brightness btwn softest quartile of Ab-injected vs control-injected
[p h] = ranksum(circMeanList{groupNums == 1}, circMeanList{groupNums == 7})

%% Make bar plots of CG intensity

paramToPlot = 'circMeanList';

figure;
p1 = bar(1, mean([circMeanList{groupNums == 7 | groupNums == 10}]), .8, 'facecolor', [0 0 .6]);
hold on;
p1e = errorbar(1, mean([circMeanList{groupNums == 7 | groupNums == 10}]), ...
    std([circMeanList{groupNums == 7 | groupNums == 10}]), 'color', 'k', 'linewidth', 2);

p2 = bar(2, mean([circMeanList{groupNums == 8 | groupNums == 11}]), .8, 'facecolor', [0 .6 0]);
hold on;
p2e = errorbar(2, mean([circMeanList{groupNums == 8 | groupNums == 11}]), ...
    std([circMeanList{groupNums == 8 | groupNums == 11}]), 'color', 'k', 'linewidth', 2);

% make width same size as p3e
h2 = get(p2e, 'children');
x2 = get(h2, 'xdata');
x2(2) = {[2 2 NaN 1.97 2.03 NaN 1.97 2.03 NaN]};
set(h2(2), 'xdata', x2{2});
set(p2e, 'children', h2);

p3 = bar(3, mean([circMeanList{groupNums == 9 | groupNums == 12}]), .8, 'facecolor', [0 0 .6]);
hold on;
p3e = errorbar(3, mean([circMeanList{groupNums == 9 | groupNums == 12}]), ...
    std([circMeanList{groupNums == 9 | groupNums == 12}]), 'color', 'k', 'linewidth', 2);

% make width same size as p3e
h1 = get(p1e, 'children');
x1 = get(h1, 'xdata');
x1(2) = {[1 1 NaN 0.97 1.03 NaN 0.97 1.03 NaN]};
set(h1(2), 'xdata', x1{2});
set(p1e, 'children', h1);

set(gca, 'fontsize', 14);
set(gca, 'xtick', [1 2 3])
xlim([.5 3.5]);
ylim([0 .15]);
set(gca, 'xticklabel', {'Too soft', 'Viable', 'Too stiff'});
ylabel('cytoplasm brightness (a.u.)');
title(sprintf('Unreleased cortical granules after fertilization'));
grid on;

%% Figure 4C in paper: CG brightness of too soft, too stiff, etc.
% Make bar plots combining CG brightness data from B6 mice (too soft or viable)
% and brightness data from microinjection controls (CBA, too soft, viable, or too stiff)

% Group1: embryos from B6 mice (n = 36)
load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Embryo\cg_imaging\embryoInfoCG.mat');
paramToPlot = 'circMeanList';


circMeanListA = embryoInfo.circMeanList;
circStdListA = embryoInfo.circStdList;
circDiffListA = embryoInfo.circDiffList;
radMeanListA = embryoInfo.radMeanList;
brightParamListA = embryoInfo.brightParamList;
viabilityA = embryoInfo.viability;

colorVec = [0 0 .6; ... % nonviable
            0 .6 0];    % viable
        
paramTooSoftA = eval([paramToPlot 'A(viabilityA == 0)']);
paramViableA = eval([paramToPlot 'A(viabilityA == 1)']);

% normalize group A params to mean viable embryo brightness
paramTooSoftA = paramTooSoftA / mean(paramViableA);
paramViableA = paramViableA / mean(paramViableA);
paramTooStiffA = []; % none that were too stiff in this data set


% Group2: CBA mice, microinjection experiment controls
fileToSave = 'C:/Users/Livia/Desktop/IVF/Processed Data/Mouse Embryo/microinjection_data/embryoInfoMicroinjection.mat';
load(fileToSave);
circMeanListB = embryoInfo.circMeanList;
circStdListB = embryoInfo.circStdList;
circDiffListB = embryoInfo.circDiffList;
radMeanListB = embryoInfo.radMeanList;
brightParamListB = embryoInfo.brightParamList;
groupNums = 1:12;

paramTooSoftB = eval(['[' paramToPlot 'B{groupNums == 7 | groupNums == 10}]']);
paramViableB = eval(['[' paramToPlot 'B{groupNums == 8 | groupNums == 11}]']);
paramTooStiffB = eval(['[' paramToPlot 'B{groupNums == 9 | groupNums == 12}]']);

% normalize group B params to mean viable embryo brightness
paramTooSoftB = paramTooSoftB / mean(paramViableB);
paramTooStiffB = paramTooStiffB / mean(paramViableB);
paramViableB = paramViableB / mean(paramViableB);

paramTooSoft = [paramTooSoftA paramTooSoftB];
paramViable = [paramViableA paramViableB];
paramTooStiff = [paramTooStiffA paramTooStiffB];


figure;
p1 = bar(1, mean(paramTooSoft), .8, 'facecolor', [0 0 .6]);
hold on;
p1e = errorbar(1, mean(paramTooSoft), std(paramTooSoft), 'color', 'k', 'linewidth', 2);

p2 = bar(2, mean(paramViable), .8, 'facecolor', [0 .6 0]);
hold on;
p2e = errorbar(2, mean(paramViable), std(paramViable), 'color', 'k', 'linewidth', 2);

% make width same size as p3e
h2 = get(p2e, 'children');
x2 = get(h2, 'xdata');
x2(2) = {[2 2 NaN 1.97 2.03 NaN 1.97 2.03 NaN]};
set(h2(2), 'xdata', x2{2});
set(p2e, 'children', h2);

p3 = bar(3, mean(paramTooStiff), .8, 'facecolor', [0 0 .6]);
hold on;
p3e = errorbar(3, mean(paramTooStiff), std(paramTooStiff), 'color', 'k', 'linewidth', 2);

% make width same size as p3e
h1 = get(p1e, 'children');
x1 = get(h1, 'xdata');
x1(2) = {[1 1 NaN 0.97 1.03 NaN 0.97 1.03 NaN]};
set(h1(2), 'xdata', x1{2});
set(p1e, 'children', h1);

set(gca, 'fontsize', 14);
set(gca, 'xtick', [1 2 3])
xlim([.5 3.5]);
ylim([0 8]);
% set(gca, 'xticklabel', {'Too soft', 'Viable', 'Too stiff'});
set(gca, 'xticklabel', {'', '', ''});
%ylabel('cytoplasm brightness (a.u.)');
title(sprintf('Unreleased cortical granules after fertilization'));
grid on;

% test for normality
[h p] = lillietest(paramTooSoft)
[h p] = lillietest(paramViable)
[h p] = lillietest(paramTooStiff)

% test for differences between nonviable and viable
% using Wilcoxon rank sum test
[p h] = ranksum(paramTooSoft, paramViable)
[p h] = ranksum(paramTooStiff, paramViable)














