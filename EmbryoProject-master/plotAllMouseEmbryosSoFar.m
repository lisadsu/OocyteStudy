% plot all mouse embryos so far
% Livia Zarnescu
% 1-27-15

clear all;
%close all;

filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse embryo\';

% don't use auto plot for 11-19 or 11-30
% whichToPlot(1) is 6-21-12 (fresh)
% whichToPlot(2) is 6-27-12 (frozen)
% whichToPlot(3) is 7-12-12 (fresh)
% whichToPlot(4) is 7-26-12 (frozen)
% whichToPlot(5) is 10-15-12 (frozen)
% whichToPlot(6) is 10-24-12 (fresh)
% whichToPlot(7) is 11-1-12 (frozen)
% whichToPlot(8) is 11-5-12 (fresh)
% whichToPlot(9) is 11-19-12 (fresh)
% whichToPlot(10) is 11-30-12 (fresh)
% whichToPlot(11) is 12-13-12 (fresh)
% whichToPlot(12) is 12-20-12 (fresh)
% whichToPlot(13) is 2-11-13 (fresh)
% whichToPlot(14) is 5-7-13 (fresh, part used for transfer) ... don't plot
% whichToPlot(15) is 6-14-13 (fresh)
% whichToPlot(16) is 7-1-13 (fresh)
% whichToPlot(17) is 7-12-13 (fresh, part used for RNA-seq) ... don't plot this one because they all lived
% whichToPlot(18) is 8-9-13 (fresh)
% whichToPlot(19) is 9-23-13 (fresh, used for RNA-seq) ... don't plot
% whichToPlot(20) is 2-19-14 (fresh, some used for transfer)
% whichToPlot(21) is 2-20-14 (fresh, some used for transfer)
% whichToPlot(22) is 6-17-14 (fresh, used for cortical granule staining)
% whichToPlot(23) is 6-27-14 (fresh, used for cortical granule staining)
% whichToPlot(24) is 7-27-14 (fresh, measured after IVF, CBA strain)
% whichToPlot(25) is 8-3-14 (fresh, measured after IVF, CBA strain)
% whichToPlot(26) is 8-14-14 (fresh, measured after IVF, CBA strain)
% whichToPlot(27) is 9-5-14 (fresh, measured after IVF, CBA strain)
% whichToPlot(28) is 9-24-14 (fresh, measured after microinj, IVF, CBA strain)
% whichToPlot(29) is 11-26-14 (fresh, measured after IVF, CBA strain)

whichToPlot = [zeros(1,8) ones(1,5) ones(1,10) 0 0 0 0 0 0];
plotAuto = [zeros(1,8) ones(1, length(whichToPlot)-8)];
pixelConv = [54*ones(1,10) 108*ones(1,length(whichToPlot)-10)];
plotCurves = 0;

% ===== MAKE SURE NEW MORPHOLOGY DATA IS SAVED =====
% =========== FILL IN HERE =========================
saveNewMorphology('mouse embryo');
load('morphologyMouseEmbryo.mat');

% morphology results of all experiments so far
% =========== FILL IN HERE ======================
morphology = cell(1,length(whichToPlot));
morphology{1} = morphology621;
morphology{2} = morphology627;
morphology{3} = morphology712;
morphology{4} = morphology726;
morphology{5} = morphology1015;
morphology{6} = morphology1024;
morphology{7} = morphology111;
morphology{8} = morphology115;
morphology{9} = morphology1119;
morphology{10} = morphology1130;
morphology{11} = morphology1213;
morphology{12} = morphology1220;
morphology{13} = morphology211;
morphology{14} = morphology5_7_13;
morphology{15} = morphology6_14_13;
morphology{16} = morphology7_1_13;
morphology{17} = morphology7_12_13;
morphology{18} = morphology8_9_13;
morphology{19} = morphology9_23_13;
morphology{20} = morphology2_19_14;
morphology{21} = morphology2_20_14;
morphology{22} = morphology6_17_14;
morphology{23} = morphology6_27_14;
morphology{24} = morphology7_27_14;
morphology{25} = morphology8_3_14;
morphology{26} = morphology8_14_14;
morphology{27} = morphology9_5_14;
morphology{28} = morphology9_24_14;
morphology{29} = morphology11_26_14;

% time lapse parameters of all experiments so far
% =========== FILL IN HERE ======================
timeLapse = cell(1,length(whichToPlot));
timeLapse{1} = timeLapse621;
timeLapse{2} = timeLapse627;
timeLapse{3} = timeLapse712;
timeLapse{4} = timeLapse726;
timeLapse{5} = timeLapse1015;
timeLapse{6} = timeLapse1024;
timeLapse{7} = timeLapse111;
timeLapse{8} = timeLapse115;
timeLapse{9} = timeLapse1119;
timeLapse{10} = timeLapse1130;
timeLapse{11} = timeLapse1213;
timeLapse{12} = timeLapse1220;
timeLapse{13} = timeLapse211;
timeLapse{14} = timeLapse5_7_13;
timeLapse{15} = timeLapse6_14_13;
timeLapse{16} = timeLapse7_1_13;
timeLapse{17} = timeLapse7_12_13;
timeLapse{18} = timeLapse8_9_13;
timeLapse{19} = timeLapse9_23_13;
timeLapse{20} = timeLapse2_19_14;
timeLapse{21} = timeLapse2_20_14;
timeLapse{22} = timeLapse6_17_14;
timeLapse{23} = timeLapse6_27_14;
timeLapse{24} = timeLapse7_27_14;
timeLapse{25} = timeLapse8_3_14;
timeLapse{26} = timeLapse8_14_14;
timeLapse{27} = timeLapse9_5_14;
timeLapse{28} = timeLapse9_24_14;
timeLapse{29} = timeLapse11_26_14;

% dates of all experiments so far
% =========== FILL IN HERE ======================
dateList = cell(1,length(whichToPlot));
dateList{1} = '6-21';
dateList{2} = '6-27';
dateList{3} = '7-12';
dateList{4} = '7-26';
dateList{5} = '10-15';
dateList{6} = '10-24';
dateList{7} = '11-1';
dateList{8} = '11-5';
dateList{9} = '11-19';
dateList{10} = '11-30';
dateList{11} = '12-13';
dateList{12} = '12-20';
dateList{13} = '2-11';
dateList{14} = '5-7-13';
dateList{15} = '6-14-13';
dateList{16} = '7-1-13';
dateList{17} = '7-12-13';
dateList{18} = '8-9-13';
dateList{19} = '9-23-13';
dateList{20} = '2-19-14';
dateList{21} = '2-20-14';
dateList{22} = '6-17-14';
dateList{23} = '6-27-14';
dateList{24} = '7-27-14';
dateList{25} = '8-3-14';
dateList{26} = '8-14-14';
dateList{27} = '9-5-14';
dateList{28} = '9-24-14';
dateList{29} = '11-26-14';

dateList2 = cell(1,length(whichToPlot));
for j = 1:length(whichToPlot)
    dateU = dateList{j};
    dateUI = strfind(dateU, '-');
    dateU(dateUI) = '_';
    dateList2{j} = dateU;
end

% for making legend to scatter plot
legendList = cell(1,3);
colorMat = [];
mList = [];
k0list = [];
k1list = [];
n0list = [];
taulist = [];
elonglist = [];
n1list = [];
gsum = [];
rsum = [];
timeLapseOut = [];
aspirationDataV = []; % keeps all aspiration curves (data, not fit)
aspirationDataN = []; % keeps all aspiration curves (data, not fit)

for jNum = 1:length(whichToPlot)
    if whichToPlot(jNum)
        jNum
        cCurr = zeros(length(morphology{jNum}), 3);
        s = size(cCurr(morphology{jNum} == 2, :));
        cCurr(morphology{jNum} == 2, :) = repmat([0 0 .6], s(1), 1);
        s = size(cCurr(morphology{jNum} == 3, :));
        cCurr(morphology{jNum} == 3, :) = repmat([0 0 .6], s(1), 1);
        s = size(cCurr(morphology{jNum} == 4, :));
        cCurr(morphology{jNum} == 4, :) = repmat([0 .6 0], s(1), 1);
        s = size(cCurr(morphology{jNum} == 1, :));
        cCurr(morphology{jNum} == 1, :) = repmat([1 0 0], s(1), 1);
        s = size(cCurr(morphology{jNum} == 5, :));
        cCurr(morphology{jNum} == 5, :) = repmat([.7 .7 0], s(1), 1);
        s = size(cCurr(morphology{jNum} == 6, :));
        cCurr(morphology{jNum} == 6, :) = repmat([.3 .3 .9], s(1), 1);
        s = size(cCurr(morphology{jNum} == 7, :));
        cCurr(morphology{jNum} == 7, :) = repmat([.3 .9 .3], s(1), 1);
        
        colorMat = [colorMat ; cCurr];
        currDateSoFar = 0;
        
        if ~isempty(timeLapse{jNum})
            timeLapseOut = [timeLapseOut; timeLapse{jNum}];
        else
            timeLapseOut = [timeLapseOut; NaN*ones(length(morphology{jNum}),3)];
        end
        
        for iNum = 1:length(morphology{jNum})
            
            if morphology{jNum}(iNum) == 3
                morphology{jNum}(iNum) = 2;
            end
            
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
                    aspirationDataN = ydata(1:37); % min length
                elseif morphology{jNum}(iNum) == 4 && isempty(legendList{2})
                    currPoint = struct('k1', k1, 'n1', n1, 'tau', tau);
                    legendList{2} = currPoint;
                    gsum = yfit;
                    aspirationDataV = ydata(1:37); % min length
                elseif morphology{jNum}(iNum) == 3 && isempty(legendList{3})
                    currPoint = struct('k1', k1, 'n1', n1, 'tau', tau);
                    legendList{3} = currPoint;
                elseif morphology{jNum}(iNum) == 2
                    rsum = [rsum; yfit];
                    aspirationDataN = [aspirationDataN; ydata(1:37)];
                elseif morphology{jNum}(iNum) == 4
                    gsum = [gsum; yfit];
                    aspirationDataV = [aspirationDataV; ydata(1:37)];
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
                
    end
end

if plotCurves
    % k0list(35) = mean(k0list([1:34 36:end]));
    xlim([0 .45]);
    % ylim([0 7.5*10^-5]);
    
    set(gca, 'FontSize', 14);
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
    title('Aspiration Depth of Mouse Embryos');
end

%% Scatter plot

% normalize variables first

numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) & ...
    ~isnan(mList);% & ~isnan(timeLapseOut(:,1))';
mN = mList(numsToPlot);
k1N = k1list(numsToPlot);
n1N = n1list(numsToPlot);
tN = taulist(numsToPlot);
k0N = k0list(numsToPlot);
tLN = timeLapseOut(numsToPlot,:);
n0N = n0list(numsToPlot);
eN = elonglist(numsToPlot);

T1 = tLN(:,1);
T2 = tLN(:,2);
T3 = tLN(:,3);
%
figure;
h = scatter3(k1N, n1N, k0N, 80, colorMat(numsToPlot,:), 'filled');
% h = scatter3(T1, T2, T3, 100, colorMat(numsToPlot,:), 'filled');

set(h, 'Marker', 'o');
set(gca, 'FontSize', 14);
title('3D scatter plot of parameters');
% title('Mouse Embryo Cell Cycle Parameters');
xlabel('k1 parameter');
ylabel('n1 parameter');
zlabel('k0 parameter');
% xlabel('cell cycle #1');
% ylabel('cell cycle #2');
% zlabel('cell cycle #3');
% axis([min(k1N) max(k1N) min(tN) max(tN)]);
% axis([min(k1N) max(k1N) min(n1N) max(n1N) min(k0N) max(k0N)]);
% axis([min(T1) max(T1) min(T2) max(T2) min(T3) max(T3)]);
% axis([0 1 0 1 0 1]);
view(-23,8);
camlight right;
drawnow;
hold on;
set(gca, 'xscale', 'linear');
set(gca, 'yscale', 'log');
set(gca, 'zscale', 'log');

mVal = 5;
%aAll = (1:length(k1list))';
aAll = (1:80)';
a = aAll(numsToPlot); %aAll(numsToPlot & mList == mVal);% & mList == mVal);% & c2);
%a = a(a > 65 & a < 85); %a(a > 84);
%a = a(a > 0 & a < 35); %a(a > 34 & a < 66);

b = num2str(a);
c = cellstr(b);
dx = -0.002; dy = 0.01; dz = .005; % displacement so the text does not overlay the data points
% text(k1N(mN == mVal & c3)+dx, n1N(mN == mVal & c3)+dy, k0N(mN == mVal & c3)+dz, c);
% text(k1list(end-104+a)+dx, n1list(end-104+a)+dy, k0list(end-104+a)+dz, c);
text(k1N(mN > mVal)+dx, n1N(mN > mVal)+dy, k0N(mN > mVal)+dz, c);
%text(k1N+dx, n1N+dy, k0N+dz, c);

hold on;

%% 
% ======== MAKE LEGEND =========
handleList = [];
labelList = cell(1,2);
colorList = [0 0 .6; ...    % dark blue
    %              0 .6 .6; ...
    0 .6 0];       % green
%
mV = [2 3 4];
for i = 1:2
    
    %     currHandle = scatter(k1N(find(mN == mV(i),1,'first')), ...
    %         n1N(find(mN == mV(i),1,'first')), ...
    %         100, colorList(i,:), 'filled');
    if ~isempty(legendList{i})
        currHandle = scatter3(k1N(find(mN == mV(i),1,'first')), ...
            n1N(find(mN == mV(i),1,'first')), k0N(find(mN == mV(i),1,'first')), ...
            80, colorList(i,:), 'filled');
        
        %     currHandle = scatter3(T1(find(mN == mV(i),1,'first')), ...
        %         T2(find(mN == mV(i),1,'first')), T3(find(mN == mV(i),1,'first')),...
        %         100, colorList(i,:), 'filled');
        
        handleList = [handleList currHandle];
        
        if i == 1
            labelList{1} = 'no blastocyst';
        elseif i == 2
            labelList{2} = 'blastocyst';
            %         elseif i == 2
            %             labelList{2} = 'poor quality blastocyst';
        end
    end
end
legend(handleList, labelList, 'Location', 'North');
legend('boxoff');
legend('boxon');
legend(handleList, labelList{1}, labelList{2}, 'Location', 'North');
% legend(handleList, labelList{1}, labelList{2}, labelList{3}, ...
%     'Location', 'North');

%% fit gaussian to each parameter

[k1_mu, k1_std] = normfit(k1N(mN == 4))
[n1_mu, n1_std] = normfit(n1N(mN == 4))
[tau_mu, tau_std] = normfit(tN(mN == 4))
[k0_mu, k0_std] = normfit(k0N(mN == 4))

% "Average" good embryo
gaussMean = [k1_mu, n1_mu, tau_mu, k0_mu]; % (k1, n1, tau, k0)
gaussStd = [k1_std, n1_std, tau_std, k0_std];

% calculate likelihood of viability vs dist from average good embryo
adjViability = mN(mN < 5);
adjViability(adjViability == 2 | adjViability == 3) = 0;
adjViability(adjViability == 4) = 1;

distList = sqrt(((k1N(mN < 5) - gaussMean(1))/gaussStd(1)).^2 + ...
    ((n1N(mN < 5) - gaussMean(2))/gaussStd(2)).^2 + ...
    ((k0N(mN < 5) - gaussMean(3))/gaussStd(3)).^2);

distFromMeanTest = sqrt(((k1N(mN == 5) - gaussMean(1))/gaussStd(1)).^2 + ...
    ((n1N(mN == 5) - gaussMean(2))/gaussStd(2)).^2 + ...
    ((k0N(mN == 5) - gaussMean(3))/gaussStd(3)).^2);

% figure, scatter(log(distList), adjViability);
%

% or ... calculate posterior probability of viability based on distFromMean
[xfit, yfit] = calcPostProb(log(distList), adjViability);
yfit = smooth(yfit,50)';

paramList = [1 1 1 1];
% fit a logistic
[bestParams, fval, exitflag] = fminsearch(@(paramList) ...
    logisticFunctionFit(xfit, yfit, paramList), [1 1 1 1])

xfit2 = xfit;
yfit2 = bestParams(1) * (1 + bestParams(2) * ...
    exp(-1*bestParams(3)*(xfit - bestParams(4)))) .^ (-1);
hold on;
plot(xfit2, yfit2);
hold on, plot(xfit, yfit);

% distance from center, normalized by standard deviation
distFromMean = sqrt(((k1N(mN == 5) - gaussMean(1))/gaussStd(1)).^2 + ...
    ((n1N(mN == 5) - gaussMean(2))/gaussStd(2)).^2 + ...
    ((k0N(mN == 5) - gaussMean(4))/gaussStd(4)).^2);

log(distFromMean)
% now find expected percentage viable from all samples with mN == 5
expectedPercentageViable = mean(interp1(xfit, yfit, log(distFromMean), 'spline', 'extrap'))



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

%% Plot of Classification Performance vs Num Params

% ======= DETECT ALL BLASTS ===========
% order: k1, tau, n1, k0
numParams = [1 2 3 4];

AUC_ROC = [.8015 .8565 .8699 .8699];
AUC_PR = [.8905 .9441 .9485 .9489];

figure, plot(numParams, AUC_ROC, '-o', 'Color', [1 0 0], 'LineWidth', 2);
hold on;
plot(numParams, AUC_PR, '-o', 'Color', [0 0 1], 'LineWidth', 2);
set(gca, 'FontSize', 14);
ylim([.7 1]);
set(gca, 'xtick',[1 2 3 4]);

xlabel('Number of Mechanical Parameters');
ylabel('Area Under Curve');
title('Classification With Mechanical Parameters');
legend('AUC_R_O_C', 'AUC_P_R');

% order: TL2, TL3, TL1
numParams = [1 2 3];

AUC_ROC = [.9206 .9551 .9328];
AUC_PR = [.9220 .9771 .9639];

figure, plot(numParams, AUC_ROC, '-o', 'Color', [1 0 0], 'LineWidth', 2);
hold on;
plot(numParams, AUC_PR, '-o', 'Color', [0 0 1], 'LineWidth', 2);
set(gca, 'FontSize', 14);
ylim([.5 1]);
set(gca, 'xtick',[1 2 3]);

xlabel('Number of Cell Cycle Parameters');
ylabel('Area Under Curve');
title('Classification With Cell Cycle Parameters');
legend('AUC_R_O_C', 'AUC_P_R');

% order: TL2, tau, TL3, k1, n1, TL1, k0
numParams = [1 2 3 4 5 6 7];

AUC_ROC = [.9206 .9518 .9603 .9559 .9497 .9406 .9369];
AUC_PR = [.9220 .9764 .9810 .9776 .9756 .9718 .9700];

figure, plot(numParams, AUC_ROC, '-o', 'Color', [1 0 0], 'LineWidth', 2);
hold on;
plot(numParams, AUC_PR, '-o', 'Color', [0 0 1], 'LineWidth', 2);
set(gca, 'FontSize', 14);
ylim([.8 1]);
set(gca, 'xtick',[1 2 3 4 5 6 7]);

xlabel('Number of Parameters');
ylabel('Area Under Curve');
title('Classification With All Parameters');
legend('AUC_R_O_C', 'AUC_P_R');


% ======= DETECT GOOD QUALITY BLASTS ===========
% order: k1, tau, n1, k0
numParams = [1 2 3 4];

AUC_ROC = [.7487 .8090 .8140 .8166];
AUC_PR = [.8458 .9065 .9014 .9017];

figure, plot(numParams, AUC_ROC, '-o', 'Color', [1 0 0], 'LineWidth', 2);
hold on;
plot(numParams, AUC_PR, '-o', 'Color', [0 0 1], 'LineWidth', 2);
set(gca, 'FontSize', 14);
ylim([.7 1]);
set(gca, 'xtick',[1 2 3 4]);

xlabel('Number of Mechanical Parameters');
ylabel('Area Under Curve');
title('Classification With Mechanical Parameters');
legend('AUC_R_O_C', 'AUC_P_R');

% order: TL2, TL3, TL1
numParams = [1 2 3];

AUC_ROC = [.8743 .9049 .8970];
AUC_PR = [.8600 .9153 .9150];

figure, plot(numParams, AUC_ROC, '-o', 'Color', [1 0 0], 'LineWidth', 2);
hold on;
plot(numParams, AUC_PR, '-o', 'Color', [0 0 1], 'LineWidth', 2);
set(gca, 'FontSize', 14);
ylim([.5 1]);
set(gca, 'xtick',[1 2 3]);

xlabel('Number of Cell Cycle Parameters');
ylabel('Area Under Curve');
title('Classification With Cell Cycle Parameters');
legend('AUC_R_O_C', 'AUC_P_R');

% order: TL2, TL3, tau, k1, TL1, n1, k0
numParams = [1 2 3 4 5 6 7];

AUC_ROC = [.8743 .9049 .9135 .9148 .9154 .9102 .9012];
AUC_PR = [.8600 .9153 .9282 .9405 .9498 .9456 .9416];

figure, plot(numParams, AUC_ROC, '-o', 'Color', [1 0 0], 'LineWidth', 2);
hold on;
plot(numParams, AUC_PR, '-o', 'Color', [0 0 1], 'LineWidth', 2);
set(gca, 'FontSize', 14);
ylim([.8 1]);
set(gca, 'xtick',[1 2 3 4 5 6 7]);

xlabel('Number of Parameters');
ylabel('Area Under Curve');
title('Classification With All Parameters');
legend('AUC_R_O_C', 'AUC_P_R');


