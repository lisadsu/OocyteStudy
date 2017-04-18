% plot all human embryos so far
% Livia Zarnescu
% 1-27-15

% clear all;
% close all;
filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Human\';

% whichToPlot(1) is 3-20-13, whichToPlot(2) is 4-2-13, both cleavage stage
% 3 is 4-17-13, cleavage stage
% 4 is 4-30-13, cleavage stage, none survived bc of incubator
% 5-13-13 is 5, 2PN stage
% 5-29-13 is 6, 2PN stage
% whichToPlot(7) is 6-28-13, 2PN stage
% whichToPlot(8) is 7-30-13, 2PN stage
% whichToPlot(9) is 11-5-13, 2PN stage
% whichToPlot(10) is 11-8-13, 2PN stage
% whichToPlot(11) is 1-14-14, 2PN stage (for RNA-seq)
% whichToPlot(12) is 2-13-14, 2PN stage (for RNA-seq)
% whichToPlot(13) is 3-19-14, 2PN stage 
% whichToPlot(14) is 6-9-14, 2PN stage (for cortical granule staining)
% whichToPlot(15) is 6-23-14, 2PN stage (for cortical granule staining)

whichToPlot = [0 0 0 0 ones(1,6) 1 1 1 0 0];
plotAuto = [0 0 0 0 ones(1,11)];
showPlot = [0 0 0 0 zeros(1,11)];

if sum(showPlot) > 0
    figure;
    hold on;
end

% ===== MAKE SURE NEW MORPHOLOGY DATA IS SAVED =====
% =========== FILL IN HERE =========================
saveNewMorphology('human');
load('morphologyHuman.mat');

% morphology results of all experiments so far
% =========== FILL IN HERE ======================
morphology = cell(1,length(whichToPlot));
morphology{1} = morphology3_20_13;
morphology{2} = morphology4_2_13;
morphology{3} = morphology4_17_13;
morphology{4} = morphology4_30_13;
morphology{5} = morphology5_13_13;
morphology{6} = morphology5_29_13;
morphology{7} = morphology6_28_13;
morphology{8} = morphology7_30_13;
morphology{9} = morphology11_5_13;
morphology{10} = morphology11_8_13;
morphology{11} = morphology1_14_14;
morphology{12} = morphology2_13_14;
morphology{13} = morphology3_19_14;
morphology{14} = morphology6_9_14;
morphology{15} = morphology6_23_14;

% time lapse parameters of all experiments so far
% =========== FILL IN HERE ======================
timeLapse = cell(1,length(whichToPlot));
timeLapse{1} = timeLapse3_20_13;
timeLapse{2} = timeLapse4_2_13;
timeLapse{3} = timeLapse4_17_13;
timeLapse{4} = timeLapse4_30_13;
timeLapse{5} = timeLapse5_13_13;
timeLapse{6} = timeLapse5_29_13;
timeLapse{7} = timeLapse6_28_13;
timeLapse{8} = timeLapse7_30_13;
timeLapse{9} = timeLapse11_5_13;
timeLapse{10} = timeLapse11_8_13;
timeLapse{11} = timeLapse1_14_14;
timeLapse{12} = timeLapse2_13_14;
timeLapse{13} = timeLapse3_19_14;
timeLapse{14} = timeLapse6_9_14;
timeLapse{15} = timeLapse6_23_14;

% dates of all experiments so far
% =========== FILL IN HERE ======================
dateList = cell(1,length(whichToPlot));
dateList{1} = '3-20-13';
dateList{2} = '4-2-13';
dateList{3} = '4-17-13';
dateList{4} = '4-30-13';
dateList{5} = '5-13-13';
dateList{6} = '5-29-13';
dateList{7} = '6-28-13';
dateList{8} = '7-30-13';
dateList{9} = '11-5-13';
dateList{10} = '11-8-13';
dateList{11} = '1-14-14';
dateList{12} = '2-13-14';
dateList{13} = '3-19-14';
dateList{14} = '6-9-14';
dateList{15} = '6-23-14';

dateList2 = cell(1,length(whichToPlot));
for j = 1:length(whichToPlot)
    dateU = dateList{j};
    dateUI = strfind(dateU, '-');
    dateU(dateUI) = '_';
    dateList2{j} = dateU;
end

% for making legend to scatter plot
legendList = cell(1,5);
colorMat = [];
k0list = [];
k1list = [];
n0list = [];
taulist = [];
elonglist = [];
n1list = [];

map = [0 0 .6; 0 .6 0];
colorIdx = [];

sumCell = cell(1,5);
mList = [];

timeLapseOut = [];

for j = 1:length(whichToPlot)
    if whichToPlot(j)
        
        morphology{j}(morphology{j} > 1 & morphology{j} < 5) = 4;
        morphology{j}(morphology{j} < 2) = 1;
        
        cCurr = zeros(length(morphology{j}), 3);
        s = size(cCurr(morphology{j} == 2, :));
        cCurr(morphology{j} == 2, :) = repmat([0 0 .6], s(1), 1);
        s = size(cCurr(morphology{j} == 3, :));
        cCurr(morphology{j} == 3, :) = repmat([0 .6 .6], s(1), 1);
        s = size(cCurr(morphology{j} == 4, :));
        cCurr(morphology{j} == 4, :) = repmat([0 .6 0], s(1), 1);
        s = size(cCurr(morphology{j} == 1, :));
        cCurr(morphology{j} == 1, :) = repmat([0 0 .6], s(1), 1);
        s = size(cCurr(morphology{j} == 5, :));
        cCurr(morphology{j} == 5, :) = repmat([.7 .7 0], s(1), 1);
        
        colorMat = [colorMat ; cCurr];
            
        if ~isempty(timeLapse{j})
            timeLapseOut = [timeLapseOut; timeLapse{j}];
        else
            timeLapseOut = [timeLapseOut; NaN*ones(length(morphology{j}),3)];
        end
        
        for i = 1:length(morphology{j})
        
        embryoString = num2str(i);
                
        % if plotAuto is on, see if an automatically measured file exists.
        % If not, just use the manually measured one.
        % if plotAuto is off, just use the manually measured one
        testPath = [filePath1, dateList{j}, ' analysis human\AutoMeasure\aspiration_data_', ...
            dateList2{j}, '_human_E', embryoString, '.mat'];
        
        if plotAuto(j) && exist(testPath)
            filePath2 = [dateList{j} ' analysis human\AutoMeasure\aspiration_data_' ...
                dateList2{j}];
        else
            filePath2 = [dateList{j} ' analysis human\aspiration_data_' ...
                dateList2{j}];
        end
            
        % if file exists, plot its data
        if exist([filePath1 filePath2 '_human_E', embryoString, '.mat']) ...
                && ~isnan(morphology{j}(i))
            
            load([filePath1 filePath2 '_human_E', embryoString, '.mat']);
            
            aspiration_depth = aspiration_depth * 40 * 10^-6 / 108; % convert from pixels to meters
            mList = [mList morphology{j}(i)];
            k0list = [k0list k0];
            k1list = [k1list k1];
            n0list = [n0list n0];
            taulist = [taulist tau];
            n1list = [n1list n1];
            elonglist = [elonglist F0/(k0 + k1)];
            
            if morphology{j}(i) == 4
                colorIdx = [colorIdx 2];
            else
                colorIdx = [colorIdx 1];
            end
            
            currColor = cCurr(i,:);
            xdata = t;
            ydata = aspiration_depth;
            
            if showPlot(j)
                plot(xdata, ydata, 'ob', 'Color', currColor);
                hold on;
                plot(xfine - min(xfine), yfit, 'Color', currColor);
            end
            
            sumCell{morphology{j}(i)} = [sumCell{morphology{j}(i)} ; yfit];
            currPoint = struct('k1', k1, 'n1', n1, 'tau', tau, 'k0', k0, ...
                'c1', timeLapse{j}(i,1), 'c2', timeLapse{j}(i,2), ...
                'c3', timeLapse{j}(i,3));
            
            % get first point of each color to add to scatterplot legend
            if morphology{j}(i) == 1 && isempty(legendList{1})
                legendList{1} = currPoint;
            elseif morphology{j}(i) == 2 && isempty(legendList{2})
                legendList{2} = currPoint;
            elseif morphology{j}(i) == 3 && isempty(legendList{3})
                legendList{3} = currPoint;
            elseif morphology{j}(i) == 4 && isempty(legendList{4})
                legendList{4} = currPoint;
            elseif morphology{j}(i) == 5 && isempty(legendList(5))
                legendList{5} = currPoint;
            end
            
        else
            mList = [mList NaN];
            k0list = [k0list NaN];
            k1list = [k1list NaN];
            n0list = [n0list NaN];
            taulist = [taulist NaN];
            n1list = [n1list NaN];
            elonglist = [elonglist NaN];
            colorIdx = [colorIdx NaN];
        end
        end
        
    end
end

if sum(showPlot) > 0
    xlim([0 .75]);
    ylim([1*10^-5 5.5*10^-5]);
    set(gca, 'FontSize', 14);
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
    title('Aspiration Depth of Human Embryos');
end

elonglist = elonglist * 10^5;


% Scatter plot

numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) ...
    & ~isnan(mList) & ~isnan(elonglist);% & ~isnan(timeLapseOut(:,1))';
mN = mList(numsToPlot);
k1N = k1list(numsToPlot);
n1N = n1list(numsToPlot);
tN = taulist(numsToPlot);
k0N = k0list(numsToPlot);
eN = elonglist(numsToPlot);
n0N = n0list(numsToPlot);

cN = colorMat(numsToPlot,:);
cI = colorIdx(numsToPlot);

tLN = timeLapseOut(numsToPlot,:);
T1 = tLN(:,1);
T2 = tLN(:,2);
T3 = tLN(:,3);


figure;
colormap(map)
h = scatter3(k1N, n1N, k0N, 200, colorIdx(numsToPlot), 'filled'); %colorMat(numsToPlot,:), 'filled');

set(h, 'Marker', 'o');
set(gca, 'FontSize', 14);
title('3D scatter plot of parameters');
xlabel('k1 parameter');
ylabel('n1 parameter');
zlabel('k0 parameter');
axis([min(k1N) max(k1N) min(n1N) max(n1N) min(k0N) max(k0N)]);
set(gca, 'yscale', 'log')
set(gca, 'zscale', 'linear');
% axis([min(k1N) max(k1N) min(n1N) max(n1N)]);

view(-23,8);
camlight right;
drawnow;
hold on;

%(54:end)
mVal = 4;
% aAll = (1:22)';%
aAll =(1:length(k1list))';
aAll = aAll(~isnan(mList));
% aAll = (1:length(k1N(end-17:end)))';%(54:end)))';
% a = aAll(numsToPlot);% & mList == mVal);
% aAll = [1:9 1:9]';
b = num2str(aAll);
c = cellstr(b);
dx = -0.001; dy = 0.015; dz = .008; % displacement so the text does not overlay the data points
% text(k1N+dx, n1N+dy, k0N+dz, c);
hold on;

% aAll(mList == mVal)'

%% patch plot

figure(2);
clf;
hold on;
theta = 0:0.1:(2*pi);
colorList = [[.85 .65 .2];[.2 .6 .9]];

for i = 1:length(k1N)
    if mN(i) < 5
        
        xPatch = cos(theta)/120 + k1N(i);
        yPatch = sin(theta)/30 + n1N(i);
        patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', colorList(cI(i),:)/2, ...
            'facecolor', colorList(cI(i),:), 'facealpha', .6, 'edgealpha', .9);

%             [xPatch, yPatch, zPatch] = sphere(10);
%             xPatch = xPatch/40 + p1(i);
%             yPatch = yPatch/20 + p2(i);
%             zPatch = zPatch/2 + p3(i);
%             currPoint = surf2patch(xPatch, yPatch, zPatch);

%             h = patch(currPoint, 'edgecolor', colorMat(i,:)/2, ...
%                 'facecolor', colorMat(i,:), 'facealpha', faceAlpha, 'edgealpha', .1);
%         
        hold on;        

    end
end

grid on;
% set(gca, 'yscale', 'log')
xlabel('k1 parameter');
ylabel('n1 parameter');

%% 6. SVM figure plot

% figure(1);
% clf;
% h = scatter(p1,p2,200,colorMat,'linewidth',4);
hold on;

axisLims = [.2 .45 .2 1.1];
embryoClassifier = fitcsvm([k1N; n1N]', cI == 2, 'ResponseName', ...
    'Viability', 'KernelFunction', 'rbf', 'KernelScale', .8, 'standardize', true);
[x1Grid,x2Grid] = meshgrid(linspace(axisLims(1),axisLims(2),100),...
    linspace(axisLims(3),axisLims(4),100));
xGrid = [x1Grid(:),x2Grid(:)];
[~,scores] = predict(embryoClassifier,xGrid);
[c, hh] = contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),-.1*[1 1],'k', 'linewidth', 2);

axis(axisLims);

% set(h, 'Marker', 'o');
% set(gca, 'FontSize', 14);
xlabel('k1 parameter');
ylabel('n1 parameter');
title('');
grid on;
axis(axisLims);



%% ======== MAKE LEGEND =========
handleList = [];
labelList = cell(1,4);
colorList = [0 0 .6; ...     % red
             0 0 .6; ...    % dark blue
             0 0 .6; ...   % light blue
             0 .6 0];       % green
for i = [1 4] 
    if ~isempty(legendList{i})
        currHandle = scatter3(legendList{i}.k1, legendList{i}.n1, ...
            legendList{i}.k0, 150, colorList(i,:), 'filled');
%         currHandle = scatter(legendList{i}.k1, legendList{i}.n1, ...
%             200, colorList(i,:), 'filled');
        
        handleList = [handleList currHandle];
        if i == 1
            labelList{1} = 'no blastocyst';
%         elseif i == 2
%             labelList{2} = 'poor blast';
%         elseif i == 3
%             labelList{3} = 'average blast';
        elseif i == 4
            labelList{4} = 'blastocyst';
        end
    else
        labelList{i} = '';
    end
end
% legend(handleList, labelList, 'Location', 'NorthWest');
% legend('boxoff');
% legend('boxon');
% legend(handleList, labelList{1}, labelList{2}, labelList{3}, labelList{4}, ...
%     'Location', 'North');

labelListShort = cell(1,2);
labelListShort{1} = labelList{1};
labelListShort{2} = labelList{4};
legend(handleList, labelListShort, 'Location', 'NorthWest');
legend('boxoff');
legend('boxon');
legend(handleList, labelList{1}, labelList{4}, ...
    'Location', 'North');

%% Plot time lapse parameters

figure;
% h = scatter3(k1N, n1N, k0N, 200, colorMat(numsToPlot,:), 'filled');
h = scatter3(T1, T2, T3, 150, colorMat(numsToPlot,:), 'filled');

set(h, 'Marker', 'o');
set(gca, 'FontSize', 14);
title('3D scatter plot of parameters');
% xlabel('k1 parameter');
% ylabel('n1 parameter');
% zlabel('k0 parameter');
xlabel('cell cycle #1');
ylabel('cell cycle #2');
zlabel('cell cycle #3');
% axis([min(k1N) max(k1N) min(n1N) max(n1N) min(k0N) max(k0N)]);
axis([min(T1) max(T1) min(T2) max(T2) min(T3) max(T3)]);
% axis([min(k1N) max(k1N) min(n1N) max(n1N)]);

view(-23,8);
camlight right;
drawnow;
hold on;

%(54:end)
mVal = 4;
% aAll = (1:length(k1list(54:end)))';
aAll = (1:length(timeLapseOut(:,1)))';%(end-20:end)))';
% a = aAll(numsToPlot);% & mList == mVal);
b = num2str(aAll);
c = cellstr(b);
dx = .1; dy = -0.5; dz = 1; % displacement so the text does not overlay the data points
% text(k1list(54:end)+dx, n1list(54:end)+dy, k0list(54:end)+dz, c);
% text(timeLapseOut(:,1)+dx, timeLapseOut(:,2)+dy, timeLapseOut(:,3)+dz, c);
hold on;
xlim([0 2]);
ylim([0 30]);

[c1_mu, c1_std] = normfit(T1(mN == 4));
[c2_mu, c2_std] = normfit(T2(mN == 4));
[c3_mu, c3_std] = normfit(T3(mN == 4));

% ======== MAKE LEGEND =========
handleList = [];
labelList = cell(1,2);
% colorList = [0 .6 .6; ...     % red
%              0 0 .6; ...    % dark blue
%              0 .6 .6; ...   % light blue
%              .8 .7 0];       % green
colorList = [0 0 .6; ...     % red
             0 0 .6; ...    % dark blue
             0 0 .6; ...   % light blue
             0 .6 0];       % green
for i = 1:4 
    if ~isempty(legendList{i})
        currHandle = scatter3(legendList{i}.c1, legendList{i}.c2, ...
            legendList{i}.c3, 150, colorList(i,:), 'filled');
%         currHandle = scatter(legendList{i}.k1, legendList{i}.n1, ...
%             200, colorList(i,:), 'filled');
        
        handleList = [handleList currHandle];
        if i == 1
            labelList{1} = 'no blastocyst';
%         elseif i == 2
%             labelList{2} = 'poor blast';
%         elseif i == 3
%             labelList{3} = 'average blast';
        elseif i == 4
            labelList{2} = 'blastocyst';
        end
    else
        %labelList{i} = '';
    end
end
legend(handleList, labelList, 'Location', 'NorthWest');
legend('boxoff');
legend('boxon');
legend(handleList, labelList{1}, labelList{2}, 'Location', 'North');
% legend(handleList, labelList{1}, labelList{2}, labelList{3}, labelList{4}, ...
%     'Location', 'North');

%% fit gaussian to each parameter

% make sure to log transform n since its distr is not normal

[k1_mu, k1_std] = normfit(k1N(mN == 4));
[n1_mu, n1_std] = normfit(log(n1N(mN == 4)));
[tau_mu, tau_std] = normfit(tN(mN == 4));
[k0_mu, k0_std] = normfit(log(k0N(mN == 4)));

gaussMean = [k1_mu, n1_mu, tau_mu, k0_mu] % (k1, n1, tau, k0)
gaussStd = [k1_std, n1_std, tau_std, k0_std]

%[k1_mu, k1_std] = normfit(k1N(mN > 2))
%[n1_mu, n1_std] = normfit(n1N(mN > 2))
%[k0_mu, k0_std] = normfit(k0N(mN > 2))

% distance from center, normalized by standard deviation
distFromMean = sqrt(((k1N(mN < 5) - gaussMean(1))/gaussStd(1)).^2 + ...
    ((log(n1N(mN < 5)) - gaussMean(2))/gaussStd(2)).^2 + ...
    ((log(k0N(mN < 5)) - gaussMean(4))/gaussStd(4)).^2)

distFromMeanTest = sqrt(((k1N(mN == 5) - gaussMean(1))/gaussStd(1)).^2 + ...
    ((log(n1N(mN == 5)) - gaussMean(2))/gaussStd(2)).^2 + ...
    ((log(k0N(mN == 5)) - gaussMean(4))/gaussStd(4)).^2)

%% Make mean curves +/- standard deviation

figure; 
set(gca, 'FontSize', 14);
plot(xfine, 10^6*mean(sumCell{1},1), 'Color', [1 0 0], 'LineWidth', 2);
hold on;
% plot(xfine, 10^6*mean(sumCell{2},1), 'Color', [0 0 .6], 'LineWidth', 2);
% plot(xfine, 10^6*mean(sumCell{3},1), 'Color', [0 .6 .6], 'LineWidth', 2);
plot(xfine, 10^6*mean(sumCell{4},1), 'Color', [0 .6 0], 'LineWidth', 2);

h = legend('No Blast', 'Good Blast');
% h = legend('No Blast', 'Poor Blast', 'Medium Blast', 'Good Blast');
set(h, 'EdgeColor', [1 1 1]);
jbfill(xfine, 10^6*(mean(sumCell{1},1)+std(sumCell{1},[],1)), 10^6*(mean(sumCell{1},1)-std(sumCell{1},[],1)), ...
    [1 0 0], 'none', [], .3);
% jbfill(xfine, 10^6*(mean(sumCell{2},1)+std(sumCell{2},[],1)), 10^6*(mean(sumCell{2},1)-std(sumCell{2},[],1)), ...
%     [0 0 .6], 'none', [], .3);
% jbfill(xfine, 10^6*(mean(sumCell{3},1)+std(sumCell{3},[],1)), 10^6*(mean(sumCell{3},1)-std(sumCell{3},[],1)), ...
%     [0 .6 .6], 'none', [], .3);
jbfill(xfine, 10^6*(mean(sumCell{4},1)+std(sumCell{4},[],1)), 10^6*(mean(sumCell{4},1)-std(sumCell{4},[],1)), ...
    [0 .6 0], 'none', [], .3);

xlim([min(xfine) max(xfine)]);
xlabel('time (seconds)');
ylabel('aspiration depth (\mum)');
title('Average Aspiration Depth Curves');


%% make animation of scatter plot rotating

angles = -140:2:-50;
angles = [angles fliplr(angles)];
elevations = 40*ones(1,length(angles));
addpath('CaptureFigVid\CaptureFigVid');

CaptureFigVid([angles; elevations]', ...
    'C:\Users\Livia\Desktop\rotatingScatterPlotHuman.avi');


%% Correlate with cell cycle parameters

A = xlsread('C:\Users\Livia\Desktop\IVF\Processed Data\Human analysis\cell cycle parameters.xlsx');

figure;
% scatter3(A(:,3)', A(:,2)', A(:,4)', 200, colorMat(~isnan(k1list),:), 'filled');
scatter3( A(:,2)',  A(:,3)', A(:,4)', 200, colorMat(~isnan(k1list(1:45)),:), 'filled');
set(gca, 'FontSize', 14);
xlabel('cell cycle 1');
ylabel('cell cycle 2');
zlabel('cell cycle 3');
title('Human embryo cell cycle parameters')

a = (1:45)';
a = a(~isnan(k1list(1:45)));
b = num2str(a);
c = cellstr(b);
dx = .05; dy = 1; dz = 1; % displacement so the text does not overlay the data points
text(A(:,2)+dx, A(:,3)+dy, A(:,4)+dz, c);
hold on;


%% Plot average curves +/- standard deviation

figure;
gavg = mean(gsum, 1);
ravg = mean(rsum, 1);

gstd = std(gsum);
rstd = std(rsum);

h1 = plot(xfine, gavg, 'Color', [0 .6 0], 'LineWidth', 2);
hold on;
% plot(xfine, gavg + gstd, 'Color', [0 .6 0]);
% plot(xfine, gavg - gstd, 'Color', [0 .6 0]);
jbfill(xfine, gavg + gstd, gavg - gstd, [0 .6 0], [0 .6 0], 0, .3);
h2 = plot(xfine, ravg, 'Color', [0 0 .6], 'LineWidth', 2);
% plot(xfine, ravg + rstd, 'r');
% plot(xfine, ravg - rstd, 'r');
jbfill(xfine, ravg + rstd, ravg - rstd, [0 0 .6], [0 0 .6], 0, .3);
set(gca, 'FontSize', 12);
title('Aspiration Curves of Fresh Embryos at 1 Cell Stage');
xlabel('time (s)');
ylabel('Aspiration Depth (pixels)');
h3 = legend([h1 h2], 'Survived to Blastocyst (n = 14)', 'Did Not Survive (n = 7)', ...
    'Location', 'SouthEast');
set(h3, 'Color', 'none');
xlim([0 max(xfine)]);

%% Make box plot of parameters

V = log(n0N(mN == 4));
N = log(n0N(mN == 1));

yV = sort(V);
yN = sort(N);
medV = median(yV);
medN = median(yN);

qUV = median(yV(yV > medV));
qLV = median(yV(yV < medV));
qUN = median(yN(yN > medN));
qLN = median(yN(yN < medN));

topV = max(yV(yV < medV + 3*(qUV - qLV)));
botV = min(yV(yV > medV - 3*(qUV - qLV)));
topN = max(yN(yN < medN + 3*(qUN - qLN)));
botN = min(yN(yN > medN - 3*(qUN - qLN)));

figure(1);
clf;
set(gca, 'fontsize', 14);
p1 = patch([.6 .6 1.4 1.4 .6], [qLN qUN qUN qLN qLN] - medV, [0 0 .6], ...
    'facealpha', .3, 'edgecolor', [0 0 .6]);
set(p1, 'linewidth', 1);
line([.6 1.4], [medN medN] - medV, 'linewidth', 2, 'color', [0 0 .6]);
line([1 1], [qUN topN] - medV, 'linewidth', 2, 'color', [0 0 .6]);
line([1 1], [qLN botN] - medV, 'linewidth', 2, 'color', [0 0 .6]);
line([.6 1.4], [topN topN] - medV, 'linewidth', 2, 'color', [0 0 .6]);
line([.6 1.4], [botN botN] - medV, 'linewidth', 2, 'color', [0 0 .6]);

hold on;
p2 = patch([1.6 1.6 2.4 2.4 1.6], [qLV qUV qUV qLV qLV] - medV, [0 .6 0], ...
    'facealpha', .3, 'edgecolor', [0 .6 0]);
set(p2, 'linewidth', 1);
line([1.6 2.4], [medV medV] - medV, 'linewidth', 2, 'color', [0 .6 0]);
line([2 2], [qUV topV] - medV, 'linewidth', 2, 'color', [0 .6 0]);
line([2 2], [qLV botV] - medV, 'linewidth', 2, 'color', [0 .6 0]);
line([1.6 2.4], [topV topV] - medV, 'linewidth', 2, 'color', [0 .6 0]);
line([1.6 2.4], [botV botV] - medV, 'linewidth', 2, 'color', [0 .6 0]);

scatter(1 + .4*(rand(1,length(N)) - .5), N - medV, 50, [0 0 .6], 'filled');
scatter(2 + .4*(rand(1,length(V)) - .5), V - medV, 50, [0 .6 0], 'filled');

line([0 3], [0 0], 'linewidth', 2, 'linestyle', '--', 'color', [.5 .5 .5]);
xlim([.3 2.7]);
set(gca, 'xtick', []);
legend([p1 p2], 'nonviable', 'viable', 'location', 'northwest');
legend('boxoff');
ylabel('\eta_0 parameter');




































