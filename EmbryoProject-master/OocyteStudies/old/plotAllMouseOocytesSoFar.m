% plot all mouse oocyte data
% Livia Zarnescu 
% 1-27-15

%%
% clear all;
% close all;
filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse oocyte\';

% whichToPlot(1) is 5-14-13 (M2)
% whichToPlot(2) is 5-30-13 (M2)
% whichToPlot(3) is 2-11-14 (M2)
% whichToPlot(4) is 2-27-14 (M2)
% whichToPlot(5) is 4-8-14 (mostly GV, some M1)
% whichToPlot(6) is 4-11-14 (mostly GV)
% whichToPlot(7) is 4-24-14
% whichToPlot(8) is 7-27-14 (M2 measured before IVF)
% whichToPlot(9) is 11-26-14 (M2 measured before IVF)

whichToPlot = [0 0 1 1 1 1 1 0 1];
plotAuto = [ones(1,length(whichToPlot))];
pixelConv = [108*ones(1,length(whichToPlot))];
plotCurves = 0;

if plotCurves
    figure;
    hold on;
end

% ===== MAKE SURE NEW MORPHOLOGY DATA IS SAVED =====
% =========== FILL IN HERE =========================
saveNewMorphology('mouse oocyte');
load('morphologyMouseOocyte.mat');

% morphology results of all experiments so far
% =========== FILL IN HERE ======================
morphology = cell(1,length(whichToPlot));
morphology{1} = morphology5_14_13;
morphology{2} = morphology5_30_13;
morphology{3} = morphology2_11_14;
morphology{4} = morphology2_27_14;
morphology{5} = morphology4_8_14;
morphology{6} = morphology4_11_14;
morphology{7} = morphology4_24_14;
morphology{8} = morphology7_27_14;
morphology{9} = morphology11_26_14;

% time lapse parameters of all experiments so far
% =========== FILL IN HERE ======================
timeLapse = cell(1,length(whichToPlot));
timeLapse{1} = timeLapse5_14_13;
timeLapse{2} = timeLapse5_30_13;
timeLapse{3} = timeLapse2_11_14;
timeLapse{4} = timeLapse2_27_14;
timeLapse{5} = timeLapse4_8_14;
timeLapse{6} = timeLapse4_11_14;
timeLapse{7} = timeLapse4_24_14;
timeLapse{8} = timeLapse7_27_14;
timeLapse{9} = timeLapse11_26_14;

% dates of all experiments so far
% =========== FILL IN HERE ======================
dateList = cell(1,length(whichToPlot));
dateList{1} = '5-14-13';
dateList{2} = '5-30-13';
dateList{3} = '2-11-14';
dateList{4} = '2-27-14';
dateList{5} = '4-8-14';
dateList{6} = '4-11-14';
dateList{7} = '4-24-14';
dateList{8} = '7-27-14';
dateList{9} = '11-26-14';

dateList2 = cell(1,length(whichToPlot));
for j = 1:length(whichToPlot)
    dateU = dateList{j};
    dateUI = strfind(dateU, '-');
    dateU(dateUI) = '_';
    dateList2{j} = dateU;
end

% for making legend to scatter plot
legendList = cell(1,4);
colorMat = [];
mList = [];
k0list = [];
k1list = [];
n0list = [];
taulist = [];
elonglist = [];
n1list = [];
timeLapseOut = [];


for jNum = 1:length(whichToPlot)
    if whichToPlot(jNum)
        jNum
        cCurr = zeros(length(morphology{jNum}), 3);
        s = size(cCurr(morphology{jNum} == 2, :));
        cCurr(morphology{jNum} == 2, :) = repmat([0 0 .6], s(1), 1);
        s = size(cCurr(morphology{jNum} == 3, :));
        cCurr(morphology{jNum} == 3, :) = repmat([0 .6 .6], s(1), 1);
        s = size(cCurr(morphology{jNum} == 4, :));
        cCurr(morphology{jNum} == 4, :) = repmat([0 .6 0], s(1), 1);
        s = size(cCurr(morphology{jNum} == 1, :));
        cCurr(morphology{jNum} == 1, :) = repmat([1 0 0], s(1), 1);
        s = size(cCurr(morphology{jNum} == 5, :));
        cCurr(morphology{jNum} == 5, :) = repmat([.7 .7 0], s(1), 1);
        s = size(cCurr(morphology{jNum} == 6, :));
        cCurr(morphology{jNum} == 6, :) = repmat([0 .6 .6], s(1), 1);
        
        colorMat = [colorMat ; cCurr];     
        currDateSoFar = 0;
        
        if ~isempty(timeLapse{jNum})
            timeLapseOut = [timeLapseOut; timeLapse{jNum}];
        else
            timeLapseOut = [timeLapseOut; NaN*ones(length(morphology{jNum}),3)];
        end

        for iNum = 1:length(morphology{jNum})
            
            embryoString = num2str(iNum);
            
            % if plotAuto is on, see if an automatically measured file exists.
            % If not, just use the manually measured one.
            % if plotAuto is off, just use the manually measured one
            testPath = [filePath1, dateList{jNum}, ' analysis\AutoMeasure\aspiration_data_', ...
                dateList2{jNum}, '_E', embryoString, '.mat'];
            
            if plotAuto(jNum) && exist(testPath)
                filePath2 = [dateList{jNum} ' analysis\AutoMeasure\aspiration_data_' ...
                    dateList2{jNum}];
            else
                filePath2 = [dateList{jNum} ' analysis\aspiration_data_' ...
                    dateList2{jNum}];
            end
        
            if exist([filePath1 filePath2 '_E', embryoString, '.mat']) && ...
                    ~isnan(morphology{jNum}(iNum))
                
                currDateSoFar = currDateSoFar + 1;
                load([filePath1 filePath2 '_E', embryoString, '.mat']);
                
                aspiration_depth = aspiration_depth * 40 * 10^-6 / pixelConv(jNum); % convert from pixels to meters
                mList = [mList morphology{jNum}(iNum)];
                k0list = [k0list k0];
                k1list = [k1list k1];
                n0list = [n0list n0];
                taulist = [taulist tau];
                n1list = [n1list n1];
                elonglist = [elonglist F0/(k0 + k1)];
                
                currColor = cCurr(iNum,:);
                xdata = t;
                ydata = aspiration_depth;
                
                if plotCurves
                    plot(xdata, 10^6*ydata, 'ob', 'Color', currColor);
                    hold on;
                    plot(xfine - min(xfine), 10^6*yfit, 'Color', currColor);
                end
                
                % get first point of each color to add to scatterplot legend
                if morphology{jNum}(iNum) == 2 && isempty(legendList{1})
                    currPoint = struct('k1', k1, 'n1', n1, 'tau', tau);
                    legendList{1} = currPoint;
                    rsum = yfit;
                elseif morphology{jNum}(iNum) == 4 && isempty(legendList{2})
                    currPoint = struct('k1', k1, 'n1', n1, 'tau', tau);
                    legendList{2} = currPoint;
                    gsum = yfit;
                elseif morphology{jNum}(iNum) == 3 && isempty(legendList{3})
                    currPoint = struct('k1', k1, 'n1', n1, 'tau', tau);
                    legendList{3} = currPoint;
                elseif morphology{jNum}(iNum) == 2
                    rsum = [rsum; yfit];
                elseif morphology{jNum}(iNum) == 4
                    gsum = [gsum; yfit];
                end
                
            else
                mList = [mList NaN];
                k0list = [k0list NaN];
                k1list = [k1list NaN];
                n0list = [n0list NaN];
                taulist = [taulist NaN];
                n1list = [n1list NaN];
                elonglist = [elonglist NaN];
            end
            
            
        end
        
        currDateSoFar

    end
end

if plotCurves
    % k0list(35) = mean(k0list([1:34 36:end]));
    xlim([0 .45]);
    % ylim([0 7.5*10^-5]);
    
    set(gca, 'FontSize', 14);
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
    title('Aspiration Depth of Mouse Oocytes');
end

% Scatter plot
% normalize variables first
numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) & ~isnan(mList);
mN = mList(numsToPlot);
k1N = k1list(numsToPlot);
n1N = n1list(numsToPlot);
tN = taulist(numsToPlot);
k0N = k0list(numsToPlot);

figure;
h = scatter(k1N, n1N, 150, colorMat(numsToPlot,:), 'filled');

set(h, 'Marker', 'o');
set(gca, 'FontSize', 16);
title('Oocyte Mechanics Reflect Maturation');
xlabel('k1 parameter');
ylabel('n1 parameter');
set(gca, 'yscale', 'log');
hold on;
grid on;

% aAll = (1:length(k1list(end-40:end)))';
aAll = (1:length(k1list))';
% a = aAll(numsToPlot);
% a = (1:20)';
b = num2str(aAll);
c = cellstr(b);
dx = -0.001; dy = 0.01; dz = .005; % displacement so the text does not overlay the data points
hold on;

% ======== MAKE LEGEND =========
handleList = [];
labelList = cell(1,3);
mLegendList = [1 2 4];
colorList = [1 0 0; ... % red
            0 0 .6; ...    % dark blue
             0 .6 0];       % green
for i = 1:3
    
    currHandle = scatter3(k1N(find(mN == mLegendList(i),1,'first')), ...
        n1N(find(mN == mLegendList(i),1,'first')), ...
        tN(find(mN == mLegendList(i),1,'first')), ...
        150, colorList(i,:), 'filled');
        handleList = [handleList currHandle];
        
        if i == 1
            labelList{1} = 'not fertilized';
        elseif i == 2
            labelList{2} = 'fertilized, no blast';
        elseif i == 3
            labelList{3} = 'fertilized, blast';
        end
end
legend(handleList, labelList, 'Location', 'North');
legend('boxoff');
legend('boxon');
% legend(handleList, labelList{1}, labelList{2}, labelList{3}, labelList{4}, ...
%     'Location', 'North');



%% Bar plot of different maturation stages (used to generate Figure 4F in paper)

figure;
p3 = bar(20, mean(k1N(mN == 2 | mN == 4)), 4, 'facecolor', [.1 .6 .6]);
hold on;
p3e = errorbar(20, mean(k1N(mN == 2 | mN == 4)), ...
    std(k1N(mN == 2 | mN == 4)), 'color', 'k', 'linewidth', 2);
setErrorBar(p3e, 20, .3);

p2 = bar(5, mean(k1N(mN == 3)), 4, 'facecolor', [.1 .6 .6]);
hold on;
p2e = errorbar(5, mean(k1N(mN == 3)), ...
    std(k1N(mN == 3)), 'color', 'k', 'linewidth', 2);
setErrorBar(p2e, 5, .3);

p1 = bar(0, mean(k1N(mN == 1)), 4, 'facecolor', [.1 .6 .6]);
hold on;
p1e = errorbar(0, mean(k1N(mN == 1)), ...
    std(k1N(mN == 1)), 'color', 'k', 'linewidth', 2);
setErrorBar(p1e, 0, .3);

plot([0 5 20], [mean(k1N(mN == 1)) mean(k1N(mN == 3)) ...
    mean(k1N(mN == 2 | mN == 4))], 'linewidth', 2, 'color', [0 0 0]);

set(gca, 'fontsize', 14);
set(gca, 'xtick', [0 5 10 15 20])
xlim([-3 23]);
ylim([0 .15]);
set(gca, 'xticklabel', {'0', '5', '10', '15', '20'});
ylabel('k_1 parameter (stiffness)');
title(sprintf('Oocyte mechanics reflect maturation'));
grid on;

% test for normality. Third one is just under threshold for normality
[h p] = lillietest(k1N(mN == 2 | mN == 4))
[h p] = lillietest(k1N(mN == 3))
[h p] = lillietest(k1N(mN == 1))

% use ranksum test to test for differencse between stages
[p h] = ranksum(k1N(mN == 2 | mN == 4), k1N(mN == 3))
[p h] = ranksum(k1N(mN == 1), k1N(mN == 3))


%% make animation of scatter plot rotating

angles = 120:2:160;
angles = [angles fliplr(angles)];
elevations = 22*ones(1,length(angles));
addpath('CaptureFigVid\CaptureFigVid');

CaptureFigVid([angles; elevations]', ...
    'C:\Users\Livia\Desktop\rotatingScatterPlotHuman3.avi');

%% Plot +/- standard deviations

figure; 
set(gca, 'FontSize', 14);
plot(xfine, 10^6*mean(gsum,1), 'Color', [0 .6 0], 'LineWidth', 2);
hold on;
plot(xfine, 10^6*mean(rsum,1), 'Color', [0 0 .6], 'LineWidth', 2);
h = legend('Blastocyst', 'No Blastocyst');
set(h, 'EdgeColor', [1 1 1]);
jbfill(xfine, 10^6*(mean(gsum,1)+std(gsum,[],1)), 10^6*(mean(gsum,1)-std(gsum,[],1)), ...
    [0 .6 0], 'none', [], .3);
jbfill(xfine, 10^6*(mean(rsum,1)+std(rsum,[],1)), 10^6*(mean(rsum,1)-std(rsum,[],1)), ...
    [0 0 .6], 'none', [], .3);

xlim([min(xfine) max(xfine)]);
xlabel('time (seconds)');
ylabel('aspiration depth (\mum)');
title('Average Aspiration Depth Curves');
