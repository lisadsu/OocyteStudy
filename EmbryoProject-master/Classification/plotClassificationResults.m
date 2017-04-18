% plot results of embryo classification
% mTest is the ground truth data in the test data set
% m_predict is the predicted embryo viability as output by the classifier
% paramsTest is a nx3 matrix with the parameters of the embryos in the test
% set
% fig_handle is just the handle to the figure to plot in
% inputMethod = 1 if ground truth is known and classifier is being
% cross-validated. InputMethod = 2 if making predictions for new embryos

function plotClassificationResults(mTest, decDist, paramsTest, fig_handle, ...
    embryoClassifier, axisLims)

m_predict = (decDist < -.5);

% Plot Results for Test data set
% green and blue for actual
ColorMat = zeros(length(mTest), 3);
num_pos = length(mTest(mTest == 1));
ColorMat(mTest == 1, :) = repmat([0 .6 0], num_pos, 1);
ColorMat(mTest == 0, :) = repmat([0 0 .6], length(mTest) - num_pos, 1);

% yellow and light blue for predicted ones
ColorPredict = zeros(length(m_predict), 3);
num_pos = length(m_predict(m_predict == 1));
ColorPredict(m_predict' == 1, :) = repmat([.7 .7 0], num_pos, 1);
ColorPredict(m_predict' == 0, :) = repmat([0 .6 .6], length(m_predict) - num_pos, 1);


if size(paramsTest,2) > 2
    
    % plot individual scatterplot points so I can make legend by color
    figure(fig_handle);
    hold on;
    h = cell(1,length(mTest));
    hLegend = cell(1,2);
    
    for i = 1:length(mTest)
        h{i} = plot3(paramsTest(i,1), paramsTest(i,2), paramsTest(i,3),...
            'marker', 'o', 'markerfacecolor', ColorMat(i,:), 'markeredgecolor', ...
            ColorMat(i,:), 'markersize', 10, 'color', 'none');
        hold on;
        if isempty(hLegend{mTest(i)+1})
            hLegend{mTest(i)+1} = h{i};
        end
    end
        
    set(gca, 'FontSize', 14);
    title('Actual Viability');
    xlim([min(paramsTest(:,1)) max(paramsTest(:,1))+.1]);
    ylim([min(paramsTest(:,2)) max(paramsTest(:,2))+.1]);
    zlim([min(paramsTest(:,3)) max(paramsTest(:,3))+.1]);
    view(152,20);
    grid on;
    axis(axisLims);
    legend([hLegend{1}, hLegend{2}], {'Nonviable', 'Viable'}, 'Location', 'North');
    
    % plot test results
    figure(fig_handle.Number+1);
    hold on;
    
    h = cell(1,length(m_predict));
    hLegend = cell(1,2);
    
    for i = 1:length(m_predict)
        h{i} = plot3(paramsTest(i,1), paramsTest(i,2), paramsTest(i,3),...
            'marker', 'o', 'markerfacecolor', ColorPredict(i,:), 'markeredgecolor', ...
            ColorPredict(i,:), 'markersize', 10, 'color', 'none');
        hold on;
        if isempty(hLegend{m_predict(i)+1})
            hLegend{m_predict(i)+1} = h{i};
        end
    end

    set(gca, 'FontSize', 14);
    title('SVM Predicted Viability');
    xlim([min(paramsTest(:,1)) max(paramsTest(:,1))+.1]);
    ylim([min(paramsTest(:,2)) max(paramsTest(:,2))+.1]);
    zlim([min(paramsTest(:,3)) max(paramsTest(:,3))+.1]);
    view(152,20);
    grid on;
    axis(axisLims);
        
    legend([hLegend{1}, hLegend{2}], {'Predicted Nonviable', ...
        'Predicted Viable'}, 'Location', 'North');

    
elseif size(paramsTest,2) == 2
    
    % plot individual scatterplot points so I can make legend by color
    figure(fig_handle);
    hold on;
    h = cell(1,length(mTest));
    hLegend = cell(1,2);
    
    for i = 1:length(mTest)
        h{i} = plot(paramsTest(i,1), paramsTest(i,2),...
            'marker', 'o', 'markerfacecolor', ColorMat(i,:), 'markeredgecolor', ...
            ColorMat(i,:), 'markersize', 10, 'color', 'none');
        hold on;
        if isempty(hLegend{mTest(i)+1})
            hLegend{mTest(i)+1} = h{i};
        end
    end
        
    set(gca, 'FontSize', 14);
    title('Actual Viability');
    xlim([min(paramsTest(:,1)) max(paramsTest(:,1))+.1]);
    ylim([min(paramsTest(:,2)) max(paramsTest(:,2))+.1]);
    grid on;
    axis(axisLims(1:4));
    legend([hLegend{1}, hLegend{2}], {'Nonviable', 'Viable'}, 'Location', 'North');
    
    % plot SVM decision boundary
    [x1Grid,x2Grid] = meshgrid(linspace(axisLims(1),axisLims(2),100),...
        linspace(axisLims(3),axisLims(4),100));
    xGrid = [x1Grid(:),x2Grid(:)];
    [~,scores] = predict(embryoClassifier,xGrid);
    contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');

    
    
    % plot test results
    figure(fig_handle.Number+1);
    hold on;
    h = cell(1,length(m_predict));
    hLegend = cell(1,2);
    
    for i = 1:length(m_predict)
        h{i} = plot(paramsTest(i,1), paramsTest(i,2),...
            'marker', 'o', 'markerfacecolor', ColorPredict(i,:), 'markeredgecolor', ...
            ColorPredict(i,:), 'markersize', 10, 'color', 'none');
        hold on;
        if isempty(hLegend{m_predict(i)+1})
            hLegend{m_predict(i)+1} = h{i};
        end
    end

    set(gca, 'FontSize', 14);
    title('SVM Predicted Viability');
    xlim([min(paramsTest(:,1)) max(paramsTest(:,1))+.1]);
    ylim([min(paramsTest(:,2)) max(paramsTest(:,2))+.1]);
    grid on;
    axis(axisLims(1:4));
        
    legend([hLegend{1}, hLegend{2}], {'Predicted Nonviable', ...
        'Predicted Viable'}, 'Location', 'North');
    
    % plot SVM decision boundary
    [x1Grid,x2Grid] = meshgrid(linspace(axisLims(1),axisLims(2),100),...
        linspace(axisLims(3),axisLims(4),100));
    xGrid = [x1Grid(:),x2Grid(:)];
    [~,scores] = predict(embryoClassifier,xGrid);
    contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
    
    
    
    
end
