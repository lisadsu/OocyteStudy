% Preliminary analysis of embryo map data

clear all;
close all;

% current folder should be Classification
addpath('..');
addpath('C:\Users\Livia\Desktop\SVM code\Bruce Code\MG_ImpactDetection');
addpath('C:\Users\Livia\Desktop\SVM code\Bruce Code\MG_ImpactDetection\Livia');
rmpath('C:\Users\Livia\Desktop\IVF\Code\Matlab Code\libsvm-3.14');


%% Load data

% specify which data set to load
method = 'mouse embryo';
nGroups = 10;
paramNumsToUse = [1 2 4];

if isequal(method, 'human embryo')
    filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Human analysis\';
    saveNewMorphology('human');
    [mOut paramsOut testIndList] = loadSVMdataHuman(1, nGroups, 0, ...
        filePath1, paramNumsToUse);
elseif isequal(method, 'mouse oocyte')
    filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse oocyte analysis\';
    saveNewMorphology('mouse oocyte');
    [mOut paramsOut testIndList] = loadSVMdataOocyte(1, nGroups, 0, ...
        filePath1, paramNumsToUse);
else % mouse embryo
    filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse embryo analysis\';
    saveNewMorphology('mouse embryo');
    [mOut paramsOut testIndList] = loadSVMdata(1, nGroups, 0, ...
        filePath1, paramNumsToUse);
end

allFeatureNames = cell(1,7);
allFeatureNames{1} = 'k1';
allFeatureNames{2} = 'n1';
allFeatureNames{3} = 'tau';
allFeatureNames{4} = 'k0';
allFeatureNames{5} = 'c1';
allFeatureNames{6} = 'c2';
allFeatureNames{7} = 'c3';

%% plot in 3D


%% Forward Feature Selection

nGroups = 20;
numRepeats = 50;
plotInput = 0;

possibleFeatures = [1:5];
featureNamesUsed = allFeatureNames(possibleFeatures);

% run forward feature selection
[featureOrderFF, ROClistFF, PRlistFF] = ...
    forwardFeatureSelectionCellCycle(mOut', ...
    paramsOut(:,possibleFeatures), ...
    nGroups, numRepeats, plotInput);

featureNamesUsed(featureOrderFF)
% featureOrderFF(featureOrderFF > 14) = featureOrderFF(featureOrderFF > 14) + 17;
featureOrderFF

save('Figures\Figures 12-10-13\Detect Good Blasts\forwardFeatureSelection_mechAndc1.mat', ...
        'featureNamesUsed', 'featureOrderFF', 'ROClistFF', 'PRlistFF');

%% Brute Force Feature Selection

% columns H-W are 1:16 (cell division frames)
% columns Y-AM are 17:31 (cell division times in hrs)
% columns AN-AP are 32:34 (cell cycle params in 2010 paper)
% column AQ is 35 (normal cell cycle params?)
% columns AR-AT are 36:38 
% (36 is # of cells before compaction)
% (37 and 38 are compaction and cavitation frames)
% columns AX-BD are 39:45 
% (39 is blast timing eval, 40 is blast ICM eval, 41 is blast TE eval)
% (42 and 43 are compaction/cavitation times)
% (44 is compaction duration, 45 is blast morphology)
% column BF is 46 (day 2 fragmentation)

deltaTimesFeatures = [];
possibleFeatures = [32:34];

inputMethod = 1;
nGroups = 20;
numRepeats = 10;
plotInput = 0;

featureNamesDT = featureNamesDeltaTimes(deltaTimesFeatures);
featureNames = allFeatureNames(possibleFeatures);
featureNamesUsed = [featureNamesDT featureNames];

% This function tries all possible combinations of inputs
% May take a long time to run
[featureOrder, bestROC, bestPR, corrZ, ROClist, PRlist] = ...
    tryFeaturesCellCycle(double(groundTruth), [deltaTimes, allParams], ...
    [deltaTimesFeatures 14+possibleFeatures], ...
    inputMethod, nGroups, numRepeats, plotInput);

save('Figures Livia\Predict Chrom Norm Only\BF_blastParamsOnly.mat', ...
        'featureNamesUsed', 'featureOrder', 'bestROC', 'bestPR');

%% classify

close all;
% fig = figure(2);
% figROC = figure(5);
groundTruth = mOut';

nGroups = 20;
crossValidate = 1;
plotInput = 0;
numRepeats = 1;

Xref = 0:.01:1;
zSVM = [0 0];

prevalence = sum(groundTruth) / length(groundTruth);

for j = 1:numRepeats
    
    % make partition
    testIndList = crossvalind('Kfold', groundTruth, nGroups)';

    % cross-validate
    [~, decDist, zSVM] = classifyExisting(testIndList, double(groundTruth)', ...
        paramsOut(:,[4 5 6 7]), fig, crossValidate, nGroups, plotInput, (j == 1), zSVM);
    
    % calculate ROC curve
    [X, Y, T, AUC] = perfcurve(double(groundTruth)', -decDist, 1);

    % get Yi out with points at a standard set of locations (at Xref)
    % otherwise finding average ROC curve isn't really possible
    Yi = interpForRoc(Xref, X, Y);
    
    if j == 1
        %         Xtotal = Xi/numRepeats; % X is 1 - specificity
        Ytotal = Yi/numRepeats; % Y is sensitivity
        Atotal = AUC;
    else
        %         Xtotal = Xtotal + Xi/numRepeats;
        Ytotal = Ytotal + Yi/numRepeats;
        Atotal = Atotal + AUC;
    end
    
%     figure(figROC);
%     hold on;
%     plot(X,Y);
%     hold on;
    
end


% figure;
% set(gca, 'FontSize', 14);
% hold on;
% plot([Xref 1], [Ytotal 1], 'Color', [0 0 1], 'LineWidth', 3);
% xlabel('1 - Specificity'); ylabel('Sensitivity')
% title('Average ROC curve');
% grid on;
% axis([0 1 0 1]);

Ztotal = Ytotal*prevalence ./ (Ytotal*prevalence + Xref*(1-prevalence)); % PPV

% figure, plot(Ytotal, Ztotal, 'Color', [0 0 1], 'LineWidth', 3);
% set(gca, 'FontSize', 14);
% xlabel('Recall (sensitivity)'); ylabel('Precision (PPV)');
% title('Average PR curve');
% grid on;


Ytotal = Ytotal(~isnan(Ztotal));
Ztotal = Ztotal(~isnan(Ztotal));

AUC_ROC = Atotal/numRepeats;
AUC_PR = trapz([0 Ytotal 1], [Ztotal(1) Ztotal Ztotal(end)]);

fprintf(['\nArea Under ROC curve is ' num2str(AUC_ROC) '\n']);
fprintf(['\nArea Under PR curve is ' num2str(AUC_PR) '\n']);

%%





figure, plot(p1, 'color', 'r', 'linewidth', 2)
% hold on; plot(p2, 'color', [0 0 .6], 'linewidth', 2)
hold on; plot(p3, 'color', [0 .6 0], 'linewidth', 2)
set(gca, 'fontsize', 14)
ylim([.6 1])
xlim([1 11])
grid on
set(gca, 'xtick', 1:11)
xlabel('Number of Params')
ylabel('AUC_R_O_C')
title('Forward Feature Selection, Params After Day 3')
legend
legend('Chromosomally Normal', 'Chrom/Blast')




        