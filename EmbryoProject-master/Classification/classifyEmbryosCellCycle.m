
% inputMethod = 1 does classification on existing ground truth data to
% evaluate goodness of classifier (cross-validation)

% inputMethod = 2 does classification of new measurements with unknown
% ground truth

% Inputs:   inputMethod:    1 does cross-validation within existing data. 
%                           2 makes predictions on new data
%           nGroups:        number of groups to use for cross-validation
%           paramNumsToUse: Vector used to specify which column numbers 
%                           in allParams to use during classification. 
%                           This gives the option to just use a few 
%                           parameters at a time for feature selection
%           groundTruth:    nx1 vector of ground truth data (0 or 1)
%           allParams:      mxn matrix where m is number of embryos
%                           measured and n is the number of features for
%                           each embryo.
%           numRepeats:     Number of times to repeat cross-validation with
%                           new random partitions each time. Should be
%                           equal to 1 when inputMethod = 2.
%           plotInput:      Output plot or no? If more than 3 features are
%                           being evaluated, only the first 3 are used to
%                           generate scatter plot
%          
% 
%

function [AUC_ROC, AUC_PR, zSVM] = classifyEmbryosCellCycle(inputMethod, nGroups, ...
    paramNumsToUse, groundTruth, allParams, numRepeats, plotInput)

close all;
% clear all;

if nargin < 4
    paramNumsToUse = 1:3;
elseif nargin < 5
    numRepeats = 10;
elseif nargin < 6
    plotInput = 0;
end

% addpath('C:\Users\Livia\Desktop\IVF\Code\Matlab Code');
% addpath('C:\Users\Livia\Desktop\IVF\Code\Matlab Code\Classification');

Xref = 0:.01:1;
zSVM = [0 0];

if inputMethod == 2
    numRepeats = 1;
end

plotInput
if plotInput
%     figROC = figure;
    fig_handle = figure;
else
%     figROC = NaN;
    fig_handle = NaN;
end

prevalence = sum(groundTruth) / length(groundTruth);

for j = 1:numRepeats
    
    %% first load data to be used for classification
    
    mOut = groundTruth;
    paramsOut = allParams(:,paramNumsToUse);
    
    if inputMethod == 1
        % Create a n-fold cross-validation to compute classification error
        testIndList = crossvalind('Kfold', groundTruth, nGroups)';
    elseif inputMethod == 2
        testIndList = ones(1,length(mOut));
        testIndList(mOut > 1) = 2;
    end
    
    
    %% Make predictions about test group
    
    nGoodSelect = 12;
    nBadSelect = 12;
    
    if inputMethod == 2
        
        % make predictions about newly measured embryos
        [worstEmbryos, bestEmbryos, m_predict, decDistTest, decDistTrain] = ...
            makeNewPredictions(testIndList, mOut, paramsOut, ...
            nGoodSelect, nBadSelect, fig_handle, [1 .7]);
        
        %         worstEmbryos = sort(worstEmbryos)
        %         bestEmbryos = sort(bestEmbryos)
        
        % find likelihood of each "test" embryo surviving
        survivalLikelihoods = findSurvivalLikelihoods(decDistTest, ...
            decDistTrain, mOut, testIndList, 2);
        
        % find best and worst embryos with this method
        [survivalSort, sInd] = sort(survivalLikelihoods, 'ascend');
        worstEmbryos2 = sInd(1:nBadSelect);
        bestEmbryos2 = sInd(end-nGoodSelect+1:end);
        
        % finally, display
        fprintf('\n Embryos Least Likely to Survive:');
        worstEmbryos2
        
        fprintf('\n Survival Likelihoods of Least Viable Embryos:');
        survivalLikelihoods(worstEmbryos2)
        
        fprintf('\n Embryos Most Likely to Survive:');
        bestEmbryos2
        
        fprintf('\n Survival Likelihoods of Most Viable Embryos:');
        survivalLikelihoods(bestEmbryos2)
        
        zSVM = [0 0];
        
        
    elseif inputMethod == 1
                
        % do cross-validation on data set that already has ground truth
        [~, decDist, zSVM] = classifyExisting(testIndList, mOut, ...
            paramsOut, fig_handle, 1, nGroups, plotInput, (j == 1), zSVM);
        
        %% Plot ROC curves
        
        % get ROC curve shape
        [X, Y, T, AUC] = perfcurve(mOut, -decDist, 1);
        
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
        
%         figure(figROC);
%         hold on;
%         plot(X,Y);
%         hold on;
        
    end
    
    
end

% calculate and plot average ROC curve if doing cross-validation
if inputMethod == 1
    
    if plotInput
        figure;
        set(gca, 'FontSize', 14);
        hold on;
        plot([Xref 1], [Ytotal 1], 'Color', [0 0 1], 'LineWidth', 3);
        xlabel('1 - Specificity'); ylabel('Sensitivity')
        title('Average ROC curve');
        grid on;
        axis([0 1 0 1]);
    end
    
    Ztotal = Ytotal*prevalence ./ (Ytotal*prevalence + Xref*(1-prevalence)); % PPV
    
    if plotInput
        figure, plot(Ytotal, Ztotal, 'Color', [0 0 1], 'LineWidth', 3);
        set(gca, 'FontSize', 14);
        xlabel('Recall (sensitivity)'); ylabel('Precision (PPV)');
        title('Average PR curve');
        grid on;
        axis([0 1 0 1]);
    end
    
    Ytotal = Ytotal(~isnan(Ztotal));
    Ztotal = Ztotal(~isnan(Ztotal));
    
    AUC_ROC = Atotal/numRepeats;
    AUC_PR = trapz([0 Ytotal 1], [Ztotal(1) Ztotal Ztotal(end)]);
    
    fprintf(['\nArea Under ROC curve is ' num2str(AUC_ROC) '\n']);
    fprintf(['\nArea Under PR curve is ' num2str(AUC_PR) '\n']);
    
else
    
    AUC_ROC = NaN;
    AUC_PR = NaN;
    
end

end
