% Forward feature selection

% Livia Zarnescu
% 4/17/14

% allParams = mxn, where m is number of data points, n is number of
%             variables

inputMethod = 1;
nGroups = 10;
numRepeats = 10;
plotInput = 0;
type = 'human embryo';

% total number of features to add
featuresToUse = [5 6 7];
totalFeatures = length(featuresToUse);

% initialize 
ROClist = zeros(1,totalFeatures);
PRlist = zeros(1,totalFeatures);
featureOrder = [];
featureUsed = zeros(1,totalFeatures);
possibleFeatures = featuresToUse; %1:totalFeatures;

% forward feature selection
for i = 1:totalFeatures
    
    % this cycle will choose the ith feature
    currROC = zeros(1,totalFeatures-i+1);
    currPR = zeros(1,totalFeatures-i+1);    
    currFeatureList = possibleFeatures(~featureUsed);
    
    % loop over possible features to add
    for j = 1:totalFeatures-i+1
                
        fprintf('\n%d features so far\n', i);
        fprintf('\nCurrent feature being tested is number %d\n', ...
            currFeatureList(j));
        
        % get performance for this specific combination of features
        [currROC(j), currPR(j), currZ(:,j)] = ...
            classifyEmbryos(inputMethod, ...
            nGroups, type, [featureOrder currFeatureList(j)], ...
            numRepeats, plotInput);
        
    end
        
    % find location of max of currROC, choose that feature and mark it as
    % chosen
    bestFeature = find(currROC == max(currROC),1);
    featureOrder = [featureOrder currFeatureList(bestFeature)];
    featureUsed(bestFeature) = 1;
    
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

% get best ROC curve vals
[~, ~, ~, ~, ~, ROCy, PRy] = ...
    classifyEmbryos(1, 10, type, ...
    featureOrder(1:find(ROClist == max(ROClist))), 10, 0);
%

if isequal(type, 'mouse embryo')
    
    if totalFeatures == 4
        featureOrder_mouse_mech = featureOrder;
        ROClist_mouse_mech = ROClist;
        PRlist_mouse_mech = PRlist;
        ROCy_mouse_mech = ROCy; % used to plot ROC curve
        PRy_mouse_mech = PRy; % used to plot PR curve
        
        if (~exist('forwardFeatureSelection.mat','file'))
            save('forwardFeatureSelection.mat', 'featureOrder_mouse_mech', ...
                'ROClist_mouse_mech', 'PRlist_mouse_mech', 'ROCy_mouse_mech', ...
                'PRy_mouse_mech');
        else
            save('forwardFeatureSelection.mat', 'featureOrder_mouse_mech', ...
                'ROClist_mouse_mech', 'PRlist_mouse_mech', 'ROCy_mouse_mech', ...
                'PRy_mouse_mech', '-append');
        end
        
    elseif totalFeatures == 7
        
        featureOrder_mouse_all = featureOrder;
        ROClist_mouse_all = ROClist;
        PRlist_mouse_all = PRlist;
        ROCy_mouse_all = ROCy; % used to plot ROC curve
        PRy_mouse_all = PRy; % used to plot PR curve
        
        if (~exist('forwardFeatureSelection.mat','file'))
            save('forwardFeatureSelection.mat', 'featureOrder_mouse_all', ...
                'ROClist_mouse_all', 'PRlist_mouse_all', 'ROCy_mouse_all', ...
                'PRy_mouse_all');
        else
            save('forwardFeatureSelection.mat', 'featureOrder_mouse_all', ...
                'ROClist_mouse_all', 'PRlist_mouse_all', 'ROCy_mouse_all', ...
                'PRy_mouse_all', '-append');
        end
    elseif totalFeatures == 3
        featureOrder_mouse_cellCycle = featureOrder;
        ROClist_mouse_cellCycle = ROClist;
        PRlist_mouse_cellCycle = PRlist;
        ROCy_mouse_cellCycle = ROCy; % used to plot ROC curve
        PRy_mouse_cellCycle = PRy; % used to plot PR curve
        
        if (~exist('forwardFeatureSelection.mat','file'))
            save('forwardFeatureSelection.mat', 'featureOrder_mouse_cellCycle', ...
                'ROClist_mouse_cellCycle', 'PRlist_mouse_cellCycle', 'ROCy_mouse_cellCycle', ...
                'PRy_mouse_cellCycle');
        else
            save('forwardFeatureSelection.mat', 'featureOrder_mouse_cellCycle', ...
                'ROClist_mouse_cellCycle', 'PRlist_mouse_cellCycle', 'ROCy_mouse_cellCycle', ...
                'PRy_mouse_cellCycle', '-append');
        end
    elseif totalFeatures == 5
        featureOrder_mouse_5feat = featureOrder;
        ROClist_mouse_5feat = ROClist;
        PRlist_mouse_5feat = PRlist;
        ROCy_mouse_5feat = ROCy; % used to plot ROC curve
        PRy_mouse_5feat = PRy; % used to plot PR curve
        
        if (~exist('forwardFeatureSelection.mat','file'))
            save('forwardFeatureSelection.mat', 'featureOrder_mouse_5feat', ...
                'ROClist_mouse_5feat', 'PRlist_mouse_5feat', 'ROCy_mouse_5feat', ...
                'PRy_mouse_5feat');
        else
            save('forwardFeatureSelection.mat', 'featureOrder_mouse_5feat', ...
                'ROClist_mouse_5feat', 'PRlist_mouse_5feat', 'ROCy_mouse_5feat', ...
                'PRy_mouse_5feat', '-append');
        end
    elseif totalFeatures == 6
        featureOrder_mouse_6feat = featureOrder;
        ROClist_mouse_6feat = ROClist;
        PRlist_mouse_6feat = PRlist;
        ROCy_mouse_6feat = ROCy; % used to plot ROC curve
        PRy_mouse_6feat = PRy; % used to plot PR curve
        
        if (~exist('forwardFeatureSelection.mat','file'))
            save('forwardFeatureSelection.mat', 'featureOrder_mouse_6feat', ...
                'ROClist_mouse_6feat', 'PRlist_mouse_6feat', 'ROCy_mouse_6feat', ...
                'PRy_mouse_6feat');
        else
            save('forwardFeatureSelection.mat', 'featureOrder_mouse_6feat', ...
                'ROClist_mouse_6feat', 'PRlist_mouse_6feat', 'ROCy_mouse_6feat', ...
                'PRy_mouse_6feat', '-append');
        end
    else
        fprintf('wrong number of features chosen ... nothing saved\n');
    end
    
elseif isequal(type, 'human embryo')
    
    if totalFeatures == 4
        featureOrder_human_mech = featureOrder;
        ROClist_human_mech = ROClist;
        PRlist_human_mech = PRlist;
        ROCy_human_mech = ROCy; % used to plot ROC curve
        PRy_human_mech = PRy; % used to plot PR curve
        
        if (~exist('forwardFeatureSelection.mat','file'))
            save('forwardFeatureSelection.mat', 'featureOrder_human_mech', ...
                'ROClist_human_mech', 'PRlist_human_mech', 'ROCy_human_mech', ...
                'PRy_human_mech');
        else
            save('forwardFeatureSelection.mat', 'featureOrder_human_mech', ...
                'ROClist_human_mech', 'PRlist_human_mech', 'ROCy_human_mech', ...
                'PRy_human_mech', '-append');
        end
    elseif totalFeatures == 7
        featureOrder_human_all = featureOrder;
        ROClist_human_all = ROClist;
        PRlist_human_all = PRlist;
        ROCy_human_all = ROCy; % used to plot ROC curve
        PRy_human_all = PRy; % used to plot PR curve
        
        if (~exist('forwardFeatureSelection.mat','file'))
            save('forwardFeatureSelection.mat', 'featureOrder_human_all', ...
                'ROClist_human_all', 'PRlist_human_all', 'ROCy_human_all', ...
                'PRy_human_all');
        else
            save('forwardFeatureSelection.mat', 'featureOrder_human_all', ...
                'ROClist_human_all', 'PRlist_human_all', 'ROCy_human_all', ...
                'PRy_human_all', '-append');
        end
    elseif totalFeatures == 3
        featureOrder_human_cellCycle = featureOrder;
        ROClist_human_cellCycle = ROClist;
        PRlist_human_cellCycle = PRlist;
        ROCy_human_cellCycle = ROCy; % used to plot ROC curve
        PRy_human_cellCycle = PRy; % used to plot PR curve
        
        if (~exist('forwardFeatureSelection.mat','file'))
            save('forwardFeatureSelection.mat', 'featureOrder_human_cellCycle', ...
                'ROClist_human_cellCycle', 'PRlist_human_cellCycle', 'ROCy_human_cellCycle', ...
                'PRy_human_cellCycle');
        else
            save('forwardFeatureSelection.mat', 'featureOrder_human_cellCycle', ...
                'ROClist_human_cellCycle', 'PRlist_human_cellCycle', 'ROCy_human_cellCycle', ...
                'PRy_human_cellCycle', '-append');
        end
    elseif totalFeatures == 5
        featureOrder_human_5feat = featureOrder;
        ROClist_human_5feat = ROClist;
        PRlist_human_5feat = PRlist;
        ROCy_human_5feat = ROCy; % used to plot ROC curve
        PRy_human_5feat = PRy; % used to plot PR curve
        
        if (~exist('forwardFeatureSelection.mat','file'))
            save('forwardFeatureSelection.mat', 'featureOrder_human_5feat', ...
                'ROClist_human_5feat', 'PRlist_human_5feat', 'ROCy_human_5feat', ...
                'PRy_human_5feat');
        else
            save('forwardFeatureSelection.mat', 'featureOrder_human_5feat', ...
                'ROClist_human_5feat', 'PRlist_human_5feat', 'ROCy_human_5feat', ...
                'PRy_human_5feat', '-append');
        end
    elseif totalFeatures == 6
        featureOrder_human_6feat = featureOrder;
        ROClist_human_6feat = ROClist;
        PRlist_human_6feat = PRlist;
        ROCy_human_6feat = ROCy; % used to plot ROC curve
        PRy_human_6feat = PRy; % used to plot PR curve
        
        if (~exist('forwardFeatureSelection.mat','file'))
            save('forwardFeatureSelection.mat', 'featureOrder_human_6feat', ...
                'ROClist_human_6feat', 'PRlist_human_6feat', 'ROCy_human_6feat', ...
                'PRy_human_6feat');
        else
            save('forwardFeatureSelection.mat', 'featureOrder_human_6feat', ...
                'ROClist_human_6feat', 'PRlist_human_6feat', 'ROCy_human_6feat', ...
                'PRy_human_6feat', '-append');
        end
    else
        fprintf('wrong number of features chosen ... nothing saved\n');
    end
end









