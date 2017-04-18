% Forward feature selection

% Livia Zarnescu
% 11/13/13

% allParams = mxn, where m is number of data points, n is number of
%             variables

function [featureOrder, ROClist, PRlist] = ...
    forwardFeatureSelectionCellCycle(groundTruth, allParams, ...
    nGroups, numRepeats, plotInput)

% total number of features to add
totalFeatures = size(allParams,2);

% initialize 
ROClist = zeros(1,totalFeatures);
PRlist = zeros(1,totalFeatures);
featureOrder = [];
featureUsed = zeros(1,totalFeatures);
possibleFeatures = 1:totalFeatures;

% forward feature selection
for i = 1:totalFeatures
    
    % this cycle will choose the ith feature
    currROC = zeros(1,totalFeatures-i+1);
    currPR = zeros(1,totalFeatures-i+1);    
    currFeatureList = possibleFeatures(~featureUsed);
    
    % loop over possible features to add
    for j = 1:totalFeatures-i+1
                
        inputMethod = 1;
        fprintf('\n%d features so far\n', i);
        fprintf('\nCurrent feature being tested is number %d\n', ...
            currFeatureList(j));
        
        % get performance for this specific combination of features
        [currROC(j), currPR(j) currZ(:,j)] = ...
            classifyEmbryosCellCycle(inputMethod, ...
            nGroups, [featureOrder currFeatureList(j)], ...
            groundTruth', allParams, numRepeats, plotInput);
        
    end
        
    % find location of max of currROC, choose that feature and mark it as
    % chosen
    bestFeature = find(currROC == max(currROC),1);
    featureOrder = [featureOrder currFeatureList(bestFeature)];
    featureUsed(currFeatureList(bestFeature)) = 1;
    
    fprintf('\n%dth best feature is number %d\n', i, ...
        currFeatureList(bestFeature));
    
    ROClist(i) = currROC(bestFeature);
    PRlist(i) = currPR(bestFeature);
    
end


figure, plot(possibleFeatures, ROClist, '-o', 'Color', [1 0 0], 'LineWidth', 2);
hold on;
plot(possibleFeatures, PRlist, '-o', 'Color', [0 0 1], 'LineWidth', 2);
set(gca, 'FontSize', 14);
if totalFeatures > 1
    xlim([1 totalFeatures]);
else
    xlim([0 2]);
end
ylim([.5 1]);
set(gca, 'xtick', possibleFeatures);

xlabel('Number of Features Used');
ylabel('Area Under Curve');
title('Classification Results');
legend('AUC_R_O_C', 'AUC_P_R', 'Location', 'SouthEast');
grid on;













