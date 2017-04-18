%
%   Inputs: 
%
%           mOut        has the ground truth data for all the embryos to be
%                       used in the training set (0 for negative group, 1 
%                       for positive group) and 2 for embryos to be predicted
%       
%           paramsOut   is the list of parameters for all the embryos (all
%                       the embryos in the training and test sets together)
%
%
%
%   Outputs: 
%
%            decDistTest    is the distance from the decision boundary for
%                           each embryo in the test set, gives a measure of
%                           the confidence of the predicted viability for
%                           each
%            decDistTrain   is the distance from the decision boundary for
%                           each embryo in the training set in order to
%                           predict likelihood of surviving for new embryos

function [decDistTest, decDistTrain] = makeNewPredictions(params, m, fig_handle, plotInput)


paramsTrain = params(m < 2,:);
paramsTest = params(m == 2, :);
mTrain = m(m < 2);
mTest = m(m == 2);

embryoClassifier = fitcsvm(paramsTrain, mTrain, 'ResponseName', 'Viability', ...
    'KernelFunction', 'rbf', 'Standardize', true);
[mTrainPredict, decDistTrainPredict] = predict(embryoClassifier, paramsTrain);
decDistTrain = decDistTrainPredict(:,2);

[mTestPredict, decDistTestPredict] = predict(embryoClassifier, paramsTest);
decDistTest = decDistTestPredict(:,2);

% plot to compare actual vs predicted (optional)
if plotInput
    
    if ~ishandle(fig_handle)
        fig_handle = figure;
    end
    
    if size(params,2) < 3
        axisLims = [min(params(:,1)) max(params(:,1)) ...
            min(params(:,2)) max(params(:,2))];
    else
        axisLims = [min(params(:,1)) max(params(:,1)) ...
            min(params(:,2)) max(params(:,2)) ...
            min(params(:,3)) max(params(:,3))];
    end
    
    % plot results
    % green and blue for actual
    ColorTrain = zeros(length(mTrain), 3);
    num_pos = length(mTrain(mTrain == 1));
    ColorTrain(mTrain == 1, :) = repmat([0 .6 0], num_pos, 1);
    ColorTrain(mTrain == 0, :) = repmat([0 0 .6], length(mTrain) - num_pos, 1);
    
    % orange and light blue for predicted ones
    ColorTest = zeros(length(mTestPredict), 3);
    num_pos = length(mTestPredict(mTestPredict == 1));
    ColorTest(mTestPredict == 1, :) = repmat([.85 .65 .2], num_pos, 1);
    ColorTest(mTestPredict == 0, :) = repmat([.2 .6 .9], length(mTestPredict) - num_pos, 1);
    
    % choose 3d or 2d plot
    if size(params,2) > 2
        
        % plot individual scatterplot points so I can make legend by color
        figure(fig_handle);
        hold on;
        hTrain = cell(1,length(mTrain));
        hTest = cell(1,length(mTest));
        hLegend = cell(1,4);
        
        for i = 1:length(mTrain)
            hTrain{i} = plot3(paramsTrain(i,1), paramsTrain(i,2), paramsTrain(i,3),...
                'marker', 'o', 'markerfacecolor', ColorTrain(i,:), 'markeredgecolor', ...
                ColorTrain(i,:), 'markersize', 10, 'color', 'none');
            if isempty(hLegend{mTrain(i)+1})
                hLegend{mTrain(i)+1} = hTrain{i};
            end
        end
        
        for i = 1:length(mTestPredict)
            hTest{i} = plot3(paramsTest(i,1), paramsTest(i,2), paramsTest(i,3),...
                'marker', 'o', 'markerfacecolor', ColorTest(i,:), 'markeredgecolor', ...
                ColorTest(i,:), 'markersize', 15, 'color', 'none');
            text(paramsTest(i,1), paramsTest(i,2), paramsTest(i,3), num2str(i));
            if isempty(hLegend{mTestPredict(i)+3})
                hLegend{mTestPredict(i)+3} = hTest{i};
            end
        end
        
        set(gca, 'FontSize', 14);
        xlim([min(params(:,1)) max(params(:,1))+.1]);
        ylim([min(params(:,2)) max(params(:,2))+.1]);
        zlim([min(params(:,3)) max(params(:,3))+.1]);
        view(152,20);
        grid on;
        axis(axisLims);
        legend([hLegend{1}, hLegend{2}, hLegend{3}, hLegend{4}], ...
            {'Nonviable', 'Viable', 'Predicted nonviable', 'Predicted viable'}, ...
            'Location', 'North');
        
    else
      
        % plot individual scatterplot points so I can make legend by color
        figure(fig_handle);
        hold on;
        hTrain = cell(1,length(mTrain));
        hTest = cell(1,length(mTest));
        hLegend = cell(1,4);
        
        for i = 1:length(mTrain)
            hTrain{i} = plot(paramsTrain(i,1), paramsTrain(i,2),...
                'marker', 'o', 'markerfacecolor', ColorTrain(i,:), 'markeredgecolor', ...
                ColorTrain(i,:), 'markersize', 10, 'color', 'none');
            if isempty(hLegend{mTrain(i)+1})
                hLegend{mTrain(i)+1} = hTrain{i};
            end
        end
        
        for i = 1:length(mTestPredict)
            hTest{i} = plot(paramsTest(i,1), paramsTest(i,2),...
                'marker', 'o', 'markerfacecolor', ColorTest(i,:), 'markeredgecolor', ...
                ColorTest(i,:), 'markersize', 15, 'color', 'none');
            text(paramsTest(i,1), paramsTest(i,2), num2str(i));
            if isempty(hLegend{mTestPredict(i)+3})
                hLegend{mTestPredict(i)+3} = hTest{i};
            end
        end
        
        set(gca, 'FontSize', 14);
        xlim([min(params(:,1)) max(params(:,1))+.1]);
        ylim([min(params(:,2)) max(params(:,2))+.1]);
        grid on;
        axis(axisLims);
        legend([hLegend{1}, hLegend{2}, hLegend{3}, hLegend{4}], ...
            {'Nonviable', 'Viable', 'Predicted nonviable', 'Predicted viable'}, ...
            'Location', 'North');
        
        % plot SVM decision boundary
        [x1Grid,x2Grid] = meshgrid(linspace(axisLims(1),axisLims(2),100),...
            linspace(axisLims(3),axisLims(4),100));
        xGrid = [x1Grid(:),x2Grid(:)];
        [~,scores] = predict(embryoClassifier,xGrid);
        scoresR = reshape(scores(:,2),size(x1Grid));
        C = contour(x1Grid,x2Grid,scoresR,0*ones(1,2),'k', 'linewidth', 2);
        
        % recalculate decision boundary as distance from 0 contour
        contourDistTest = zeros(1,length(decDistTest));
        contourDistTrain = zeros(1,length(decDistTrain));
        
        for i = 1:length(decDistTest)
            paramDiffZeroContour = C' - repmat(paramsTest(i,:), size(C,2), 1);
            distFromZeroContour = sqrt(paramDiffZeroContour(:,1).^2 + paramDiffZeroContour(:,2).^2);
            contourDistTest(i) = sign(decDistTest(i))*min(distFromZeroContour);
        end
        
        for i = 1:length(decDistTrain)
            paramDiffZeroContour = C' - repmat(paramsTrain(i,:), size(C,2), 1);
            distFromZeroContour = sqrt(paramDiffZeroContour(:,1).^2 + paramDiffZeroContour(:,2).^2);
            contourDistTrain(i) = sign(decDistTrain(i))*min(distFromZeroContour);
        end
        
        % for 2D only
        decDistTest = contourDistTest';
        decDistTrain = contourDistTrain';
        
    end
    
    
end
