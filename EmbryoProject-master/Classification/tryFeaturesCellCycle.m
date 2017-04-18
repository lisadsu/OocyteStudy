% Script to try different combinations of features
% currently finds outcome of all combinations of input features
% i.e. this is the brute force approach
% may take excessive computational time

% Livia Zarnescu
% 10/2/13

% allParams = mxn, where m is number of data points, n is number of
%             variables
% possibleFeatures = which features in allParams to evaluate all possible
%                    combinations of
% inputMethod should = 1
% 

function [featureOrder, bestROC, bestPR, corrZ, ROClist, PRlist] = ...
    tryFeaturesCellCycle(groundTruth, allParams, possibleFeatures, ...
    inputMethod, nGroups, numRepeats, plotInput)

numTotalFeatures = length(possibleFeatures);
allNums = 1:(2^numTotalFeatures - 1);
totalCombinations = length(allNums)

allNumsBin = dec2bin(allNums,numTotalFeatures);

% sort by number of features used
[sumSort, indSort] = sort(sum(allNumsBin,2), 'ascend');

% This is the combination vector
allNumsBinSort = allNumsBin(indSort,:);

sumList = min(sumSort):max(sumSort);

% allCombNFeat{i} contains a logical vector with the indices of feature
% combinations with exactly i features at once
allCombNFeat = cell(1,numTotalFeatures);

for i = 1:numTotalFeatures
   allCombNFeat{i} =  (sumSort == sumList(i));
end

% 2. Evaluate classification performance on varying combinations of features

ROClist = zeros(1,totalCombinations);
PRlist = zeros(1,totalCombinations);
zList = zeros(2,totalCombinations);

for i = 1:totalCombinations
    
    i
    % get performance for this specific combination of features
    [ROClist(i), PRlist(i) zList(:,i)] = classifyEmbryosCellCycle(inputMethod, ...
        nGroups, possibleFeatures(allNumsBinSort(i,:) == '1'), groundTruth', ...
        allParams, numRepeats, plotInput);
    
end


% 3. Print out appropriate feature order and corresponding AUC's

featureOrder = cell(1,numTotalFeatures);
bestROC = zeros(1,numTotalFeatures);
bestPR = zeros(1,numTotalFeatures);
corrZ = zeros(2,numTotalFeatures);

for i = 1:numTotalFeatures
   
    % set all values with diff # of features to 0
    currROC = ROClist;
    currPR = PRlist;
    currZ = zList;
    currROC(~allCombNFeat{i}) = 0;
    currPR(~allCombNFeat{i}) = 0;
    currZ(:,~allCombNFeat{i}) = 0;
        
    % find index of highest AUC_ROC in remaining values
    bestInd = find(currROC == max(currROC), 1, 'first');
    bestROC(i) = currROC(bestInd);
    bestPR(i) = currPR(bestInd);
    corrZ(:,i) = currZ(:,bestInd);
    featuresUsed = possibleFeatures(allNumsBinSort(bestInd,:) == '1')
    
    featureOrder{i} = featuresUsed;
    
end

bestROC
bestPR

figure, plot(1:numTotalFeatures, bestROC, '-o', 'Color', [1 0 0], 'LineWidth', 2);
hold on;
plot(1:numTotalFeatures, bestPR, '-o', 'Color', [0 0 1], 'LineWidth', 2);
set(gca, 'FontSize', 14);
xlim([1 numTotalFeatures]);
ylim([.5 1]);
set(gca, 'xtick',1:numTotalFeatures);

xlabel('Number of Parameters');
ylabel('Area Under Curve');
title('Classification Results');
legend('AUC_R_O_C', 'AUC_P_R', 'Location', 'SouthEast');
grid on;

save('forwardFeatureSelection.mat', 'featureOrder', 'bestROC', 'bestPR', 'corrZ')











