% Script to try different combinations of features
% Can be used for forward feature selection

% Livia Zarnescu
% 9/12/13

% 1. populate combination vector

possibleFeatures = [1 2 3 4];
numTotalFeatures = length(possibleFeatures);
allNums = 1:(2^numTotalFeatures - 1);
totalCombinations = length(allNums);

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
ROCstd = zeros(1,totalCombinations);
PRstd = zeros(1,totalCombinations);
zList = zeros(2,totalCombinations);

for i = 1:totalCombinations
    
    i
    [ROClist(i), PRlist(i), zList(:,i), ROCstd(i), PRstd(i)] = ...
        classifyEmbryos(1,10,'mouse embryo', ...
       possibleFeatures(allNumsBinSort(i,:) == '1'), 10, 0);
    
end


% 3. Print out appropriate feature order and corresponding AUC's

featureOrder = cell(1,numTotalFeatures);
bestROC = zeros(1,numTotalFeatures);
bestROCstd = zeros(1,numTotalFeatures);
bestPR = zeros(1,numTotalFeatures);
bestPRstd = zeros(1,numTotalFeatures);
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
    bestROCstd(i) = ROCstd(bestInd);
    bestPR(i) = currPR(bestInd);
    bestPRstd(i) = PRstd(bestInd);
    corrZ(:,i) = currZ(:,bestInd);
    featuresUsed = possibleFeatures(allNumsBinSort(bestInd,:) == '1')
    
    featureOrder{i} = featuresUsed;
    
end

bestROC
bestROCstd
bestPR
bestPRstd

figure, plot(1:numTotalFeatures, bestROC, '-o', 'Color', [1 0 0], 'LineWidth', 2);
hold on;
plot(1:numTotalFeatures, bestPR, '-o', 'Color', [0 0 1], 'LineWidth', 2);
set(gca, 'FontSize', 14);
xlim([1 numTotalFeatures]);
ylim([.8 1]);
set(gca, 'xtick',1:numTotalFeatures);

xlabel('Number of Parameters');
ylabel('Area Under Curve');
title('Classification With Mechanical Parameters');
legend('AUC_R_O_C', 'AUC_P_R', 'Location', 'SouthEast');
grid on;






