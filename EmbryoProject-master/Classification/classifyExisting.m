% Perform cross-validation on existing data
%
%

function [predictOut, decDistOut, fig_handle] = classifyExisting(params, m, fig_handle, plotInput)

% ======================================================
% Perform cross-validation between the n groups
% ======================================================

% train SVM classifier with 10-fold cross-validation
% make predictions and return decision boundary distances

bestSigma = .7; % findBestSVMParams(params, m);

embryoClassifier = fitcsvm(params, m, 'ResponseName', 'Viability', ...
    'KernelFunction', 'rbf', 'Standardize', true, 'CrossVal', 'on');%, 'KernelScale', bestSigma);
[predictOut, decDist] = kfoldPredict(embryoClassifier);
decDistOut = decDist(:,2);


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
    ColorMat = zeros(length(m), 3);
    num_pos = length(m(m == 1));
    ColorMat(m == 1, :) = repmat([0 .6 0], num_pos, 1);
    ColorMat(m == 0, :) = repmat([0 0 .6], length(m) - num_pos, 1);
    
    % orange and light blue for predicted ones
    ColorPredict = zeros(length(predictOut), 3);
    num_pos = length(predictOut(predictOut == 1));
    ColorPredict(predictOut' == 1, :) = repmat([.85 .65 .2], num_pos, 1);
    ColorPredict(predictOut' == 0, :) = repmat([.2 .6 .9], length(predictOut) - num_pos, 1);
    
    % choose 3d or 2d plot
    if size(params,2) > 2
        
        % plot individual scatterplot points so I can make legend by color
        figure(fig_handle);
        subplot(1,2,1);
        hold on;
        h = cell(1,length(m));
        hLegend = cell(1,2);
        
        for i = 1:length(m)
            h{i} = plot3(params(i,1), params(i,2), params(i,3),...
                'marker', 'o', 'markerfacecolor', ColorMat(i,:), 'markeredgecolor', ...
                ColorMat(i,:), 'markersize', 10, 'color', 'none');
            hold on;
            if isempty(hLegend{m(i)+1})
                hLegend{m(i)+1} = h{i};
            end
        end
        
        set(gca, 'FontSize', 14);
        title('Actual Viability');
        xlim([min(params(:,1)) max(params(:,1))+.1]);
        ylim([min(params(:,2)) max(params(:,2))+.1]);
        zlim([min(params(:,3)) max(params(:,3))+.1]);
        view(152,20);
        grid on;
        axis(axisLims);
        legend([hLegend{1}, hLegend{2}], {'Nonviable', 'Viable'}, 'Location', 'North');
        
        
        subplot(1,2,2);
        hold on;
        
        h = cell(1,length(predictOut));
        hLegend = cell(1,2);
        
        for i = 1:length(predictOut)
            h{i} = plot3(params(i,1), params(i,2), params(i,3),...
                'marker', 'o', 'markerfacecolor', ColorPredict(i,:), 'markeredgecolor', ...
                ColorPredict(i,:), 'markersize', 10, 'color', 'none');
            hold on;
            if isempty(hLegend{predictOut(i)+1})
                hLegend{predictOut(i)+1} = h{i};
            end
        end
        
        set(gca, 'FontSize', 14);
        title('SVM Predicted Viability');
        xlim([min(params(:,1)) max(params(:,1))+.1]);
        ylim([min(params(:,2)) max(params(:,2))+.1]);
        zlim([min(params(:,3)) max(params(:,3))+.1]);
        view(152,20);
        grid on;
        axis(axisLims);
        
        legend([hLegend{1}, hLegend{2}], {'Predicted Nonviable', ...
            'Predicted Viable'}, 'Location', 'North');
        
    else
      
        % plot individual scatterplot points so I can make legend by color
        figure(fig_handle);
        subplot(1,2,1);
        hold on;
        h = cell(1,length(m));
        hLegend = cell(1,2);
        
        for i = 1:length(m)
            h{i} = plot(params(i,1), params(i,2),...
                'marker', 'o', 'markerfacecolor', ColorMat(i,:), 'markeredgecolor', ...
                ColorMat(i,:), 'markersize', 10, 'color', 'none');
            hold on;
            if isempty(hLegend{m(i)+1})
                hLegend{m(i)+1} = h{i};
            end
        end
        
        set(gca, 'FontSize', 14);
        title('Actual Viability');
        xlim([min(params(:,1)) max(params(:,1))+.1]);
        ylim([min(params(:,2)) max(params(:,2))+.1]);
        grid on;
        axis(axisLims);
        legend([hLegend{1}, hLegend{2}], {'Nonviable', 'Viable'}, 'Location', 'North');
        
        % plot SVM decision boundary
        embryoClassifierTrainAll = fitcsvm(params, m, 'ResponseName', ...
            'Viability', 'KernelFunction', 'rbf', 'Standardize', true);%, 'KernelScale', bestSigma);
        [x1Grid,x2Grid] = meshgrid(linspace(axisLims(1),axisLims(2),100),...
            linspace(axisLims(3),axisLims(4),100));
        xGrid = [x1Grid(:),x2Grid(:)];
        [~,scores] = predict(embryoClassifierTrainAll,xGrid);
        contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),-.0*ones(1,2),'k');

        
        subplot(1,2,2);
        hold on;
        
        h = cell(1,length(predictOut));
        hLegend = cell(1,2);
        
        for i = 1:length(predictOut)
            h{i} = plot(params(i,1), params(i,2),...
                'marker', 'o', 'markerfacecolor', ColorPredict(i,:), 'markeredgecolor', ...
                ColorPredict(i,:), 'markersize', 10, 'color', 'none');
            hold on;
            if isempty(hLegend{predictOut(i)+1})
                hLegend{predictOut(i)+1} = h{i};
            end
        end
        
        set(gca, 'FontSize', 14);
        title('SVM Predicted Viability');
        xlim([min(params(:,1)) max(params(:,1))+.1]);
        ylim([min(params(:,2)) max(params(:,2))+.1]);
        grid on;
        axis(axisLims);
        
        legend([hLegend{1}, hLegend{2}], {'Predicted Nonviable', ...
            'Predicted Viable'}, 'Location', 'North');
        
%         % plot SVM decision boundary
%         embryoClassifierTrainAll = fitcsvm(params, m, 'ResponseName', ...
%             'Viability', 'KernelFunction', 'rbf', 'KernelScale', bestSigma);
%         [x1Grid,x2Grid] = meshgrid(linspace(axisLims(1),axisLims(2),100),...
%             linspace(axisLims(3),axisLims(4),100));
%         xGrid = [x1Grid(:),x2Grid(:)];
%         [~,scores] = predict(embryoClassifierTrainAll,xGrid);
%         contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');

        
    end
    
    
end




