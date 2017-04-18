% analyze oocyte mechanics vs fertilization
% Livia Zarnescu Yanez
% 4-12-16


params = loadDataOocytes();
experimentType = params.experimentType;
numExperiments = params.numExperiments;
dateList = params.dateList;
numOocytes = params.numOocytes;
oocyteNums = params.oocyteNums;
fertInfo = params.fertInfo;
blastForm = params.blastForm;
hatchInfo = params.hatchInfo;
maturationEnv = params.maturationEnv;
k1ScaleFactor = params.k1ScaleFactor;
morphologyInfo = params.morphologyInfo;
measHour = params.measHour;


procDataPath = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\';

% init param vectors
fertList = [];
blastList = [];
hatchList = [];


k0list = [];
k1list = [];
n0list = [];
n1list = [];
taulist = [];
colorMat = [];
aList = [];
matEnv = [];
legendParam = [];
oocyteM = [];
measHourList = [];
eType = [];
currE = 0;
colorList = [[.2 .2 .6];[.8 .2 .8];[.6 .6 .6]];
%colorList = [[.5 .5 1];[.3 .3 .7];[.1 .1 .4]];

exptsToPlot = ones(1,12);
aText = [];
cText = {};

% load in params by participant and embryo num
for i = 1:numExperiments
    if exptsToPlot(i)
        
        currDate = dateList{i}
        currDateU = strrep(currDate,'-','_');
        
        for j = 1:numOocytes(i)
            
            currE = currE + 1;
            currDataPath = [procDataPath currDate ' analysis\AutoMeasure\aspiration_data_' ...
                currDateU '_E' num2str(oocyteNums{i}(j)) params.fileNameApp{i}{j} '.mat'];
            
            if ~exist(currDataPath, 'file')
                currDataPath = [procDataPath currDate ' analysis\aspiration_data_' ...
                    currDateU '_E' num2str(oocyteNums{i}(j)) params.fileNameApp{i}{j} '.mat'];
            end
            
            % save embryo params and color
            if exist(currDataPath, 'file') && (morphologyInfo{i}(j) > -1)% && ~isnan(fertInfo{i}(j))
                
                load(currDataPath);
                
                fertList = [fertList fertInfo{i}(j)];
                blastList = [blastList blastForm{i}(j)];
                hatchList = [hatchList hatchInfo{i}(j)];
                k0list = [k0list k0];
                k1list = [k1list (k1+k1ScaleFactor(i))];
                n0list = [n0list n0];
                n1list = [n1list (n1+k1ScaleFactor(i)*.4059/.0314)];
                taulist = [taulist tau];
                oocyteM = [oocyteM morphologyInfo{i}(j)];
                matEnv = [matEnv maturationEnv{i}(j)];
                aText = [aText currE];
                measHourList = [measHourList measHour{i}(j)];
                eType = [eType experimentType(i)];
                cText = {cText{:}, [currDate, '\_E', num2str(oocyteNums{i}(j)) params.fileNameApp{i}{j}]};
                
                aPad = padarray(A,[0,65-length(A)], 'replicate', 'post');
                aList = [aList; aPad(1:65)];
                
                %                 % color code based on fertilization and blast formation
                %                 legendParamNums = 3;
                %                 if (fertList(end) == 0)
                %                     legendParam = [legendParam 1];
                %                     colorMat = [colorMat; [.85 .65 .2]]; % no fert
                %                 elseif (hatchList(end) == 0)
                %                     legendParam = [legendParam 2];
                %                     colorMat = [colorMat; [.9 .2 .8]]; % fert, no blast
                %                 else
                %                     legendParam = [legendParam 3];
                %                     colorMat = [colorMat; [0 .2 .8]]; % fert, blast
                %                 end
                
                
                % colorcode based on fertilization time
                legendParamNums = 4;
                if measHour{i}(j) < 11
                    legendParam = [legendParam 1];
                    colorMat = [colorMat; colorList(1,:)];%[.85 .65 .2]]; % 6 hrs
                    %                 elseif measHour{i}(j) < 11
                    %                     legendParam = [legendParam 2];
                    %                     colorMat = [colorMat;  colorList(2,:)]; % 10 hrs
                elseif measHour{i}(j) < 17
                    legendParam = [legendParam 2];
                    colorMat = [colorMat;  colorList(2,:)]; % 14 hrs
                elseif measHour{i}(j) > 17
                    legendParam = [legendParam 3];
                    colorMat = [colorMat;  colorList(3,:)]; % 20 hrs
                end
                
                
                % colorcode based on maturation environment
                %                 legendParamNums = 3;
                %                 if (maturationEnv{i}(j) == 0)
                %                     legendParam = [legendParam 1];
                %                     colorMat = [colorMat; [.2 .6 .9]]; % in vivo
                %                 elseif (maturationEnv{i}(j) == 1)
                %                     legendParam = [legendParam 3];
                %                     colorMat = [colorMat; [.85 .65 .2]]; % KSOM
                %                 else
                %                     legendParam = [legendParam 2];
                %                     colorMat = [colorMat; [.9 .2 .8]]; % MM
                %                 end
                
                %                 % colorcode based on oocyte morphology
                %                 legendParamNums = 3;
                %                 if (morphologyInfo{i}(j) == 2)
                %                     legendParam = [legendParam 1];
                %                     colorMat = [colorMat; [.8 .2 .8]]; % MII
                %                 elseif (morphologyInfo{i}(j) == 1)
                %                     legendParam = [legendParam 2];
                %                     colorMat = [colorMat; [.5 .5 .5]]; % MI
                %                 else
                %                     legendParam = [legendParam 3];
                %                     colorMat = [colorMat; [0 .2 .8]]; % GV
                %                 end
                
            else
                colorMat = [colorMat; [NaN NaN NaN]];
                fertList = [fertList NaN];
                blastList = [blastList NaN];
                hatchList = [hatchList NaN];
                k0list = [k0list NaN];
                k1list = [k1list NaN];
                n0list = [n0list NaN];
                n1list = [n1list NaN];
                taulist = [taulist NaN];
                oocyteM = [oocyteM NaN];
                matEnv = [matEnv NaN];
                eType = [eType NaN];
                legendParam = [legendParam NaN];
                aList = [aList; NaN*zeros(1,65)];
                aText = [aText NaN];
                measHourList = [measHourList NaN];
                cText = {cText{:}, ''};
            end
        end
    end
end

aList = 10^6*aList;
p1 = k1list;
p2 = n1list;
p3 = oocyteM;% + rand(1,length(oocyteM))/2-.25;
%%
figure(1);
clf;
hold on;
h = cell(1,length(k1list));
hLegend = cell(1,legendParamNums);
theta = 0:0.1:(2*pi);

for i = 1:length(k1list)
    i
    if ~isnan(legendParam(i))
        i
        %         h{i} = plot(p1(i), p2(i), 'marker', 'o', 'markeredgecolor', ...
        %             colorMat(i,:), 'markersize', 12, 'color', 'none', ...
        %             'linewidth', 4);
        %
        xPatch = cos(theta)*(max(p1) - min(p1))/50 + p1(i);
        yPatch = sin(theta)*(max(p2) - min(p2))/50 + p2(i);
        
        h{i} = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', colorMat(i,:)/2, ...
            'facecolor', colorMat(i,:), 'facealpha', .5, 'edgealpha', .5);
        
        
        %         [xPatch, yPatch, zPatch] = sphere(10);
        %         xPatch = xPatch*.003 + p1(i);
        %         yPatch = yPatch*.3 + p2(i);
        %         zPatch = zPatch*.07 + p3(i);
        %         currPoint = surf2patch(xPatch, yPatch, zPatch);
        %
        %         h{i} = patch(currPoint, 'edgecolor', colorMat(i,:)/2, ...
        %             'facecolor', colorMat(i,:), 'facealpha', .3, 'edgealpha', 0);
        %
        % %         h{i} = plot3(p1(i), p2(i), p3(i), 'marker', 'o', 'markeredgecolor', ...
        % %             colorMat(i,:), 'markersize', 12, 'color', colorMat(i,:), ...
        % %             'linewidth', 4);
        
        hold on;
        if isempty(hLegend{legendParam(i)})
            i
            hLegend{legendParam(i)} = h{i};
        end
    end
end

hourList = [6 14 20];
% for i = 1:length(hourList)
%     p1Curr = mean(p1(measHourList == hourList(i)));
%     p2Curr = mean(p2(measHourList == hourList(i)));
%     p1Std = std(p1(measHourList == hourList(i)));
%     p2Std = std(p2(measHourList == hourList(i)));
%
%     hh = herrorbar(p1Curr, p2Curr, p1Std);
%     set(hh, 'color', [0 0 0], 'linewidth', 1);
%
%     hv = terrorbar(p1Curr, p2Curr, p2Std, 0);
%     set(hv, 'color', [0 0 0], 'linewidth', 1);
% end


% hLegend = {};
% for i = 1:length(hourList)
%
%     p1Curr = median(p1(measHourList == hourList(i)));
%     p2Curr = median(p2(measHourList == hourList(i)));
%     xW = .003;
%     yW = .35;
%
%     xPatch = [p1Curr - xW, p1Curr - xW, p1Curr + xW, p1Curr + xW];
%     yPatch = [p2Curr - yW, p2Curr + yW, p2Curr + yW, p2Curr - yW];
%
%     hLegend{i} = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', colorList(i,:)/2, ...
%             'facecolor', 1/4 + 3/4*colorList(i,:), 'facealpha', 1, 'edgealpha', 1);
%
% end



set(gca, 'FontSize', 14);
title('Oocyte outcomes after IVF');
grid on;
axis([.06 .155 -.5 9]);
xlabel('k_1 parameter');
ylabel('\eta_1 parameter');
zlabel('morphology');

% set(gca, 'ztick', [1 2]);
% set(gca, 'zticklabel', {'MI', 'MII'});
legend([hLegend{1}, hLegend{2}, hLegend{3}], ...
    {'No fert', 'Fert', 'Hatched blast'}, 'Location', 'NorthEast');
% legend([hLegend{1}, hLegend{2}, hLegend{3}], ...
%     {'8 hr', '14 hr', '20 hr'}, 'Location', 'NorthEast');
% legend([hLegend{1}, hLegend{2}, hLegend{3}], ...
%     {'In vivo', 'MM', 'KSOM'}, 'Location', 'NorthEast');
% legend([hLegend{1}, hLegend{2}, hLegend{3}], ...
%     {'MII', 'MI', 'GV'}, 'Location', 'NorthEast');

b = num2str(aText');
c = cellstr(b);
dx = -0.01; dy = 0.015; dz = .5; % displacement so the text does not overlay the data points
hold on;
% view(-33, 24);
% text(p1+dx, p2+dy, cText'); %, p3+dz, c);
% view(.5,90)
% view(.5,0)



%% Plot correlation between morphology and mechanical params

figure(3);
clf;

k10 = bar(1, mean(k1list(oocyteM == 0 & matEnv == 0)), .8, 'facecolor', [.2 .6 .6]);
hold on;
ek10 = terrorbar(1, mean(k1list(oocyteM == 0 & matEnv == 0)), std(k1list(oocyteM == 0 & matEnv == 0)), .1);
set(ek10, 'color', 'k', 'linewidth', 2);

k11 = bar(2, mean(k1list(oocyteM == 1 & matEnv == 0)), .8, 'facecolor', [.4 .8 .8]);
ek11 = terrorbar(2, mean(k1list(oocyteM == 1 & matEnv == 0)),std(k1list(oocyteM == 1 & matEnv == 0)), .1);
set(ek11, 'color', 'k', 'linewidth', 2);

k12 = bar(3, mean(k1list(oocyteM == 2 & matEnv == 0)), .8, 'facecolor', [.8 1 1]);
ek12 = terrorbar(3, mean(k1list(oocyteM == 2 & matEnv == 0)),std(k1list(oocyteM == 2 & matEnv == 0)), .1);
set(ek12, 'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [1 2 3])
set(gca, 'fontsize', 14);
ylim([0 .3]);
set(gca, 'xticklabel', {'GV', 'MI', 'MII'});
ylabel('k_1 parameter');
title('Effect of nuclear maturation on mechanics');
xlim([0.5 3.5]);

[p h] = ranksum(k1list(oocyteM == 2 & matEnv == 0), k1list(oocyteM == 1 & matEnv == 0))
[p h] = ranksum(k1list(oocyteM == 2 & matEnv == 0), k1list(oocyteM == 0 & matEnv == 0))


%% Plot bar chart with all combinations of maturation environment and morphology


figure(4);
clf;

% In Vivo matured
k00 = bar(0, mean(k1list(oocyteM == 0 & matEnv == 0)), .25, 'facecolor', [.2 .6 .6]);
hold on;
ek00 = terrorbar(0, mean(k1list(oocyteM == 0 & matEnv == 0)), std(k1list(oocyteM == 0 & matEnv == 0)), .1);
set(ek00, 'color', 'k', 'linewidth', 2);

k10 = bar(.25, mean(k1list(oocyteM == 1 & matEnv == 0)), .25, 'facecolor', [.4 .8 .8]);
hold on;
ek10 = terrorbar(.25, mean(k1list(oocyteM == 1 & matEnv == 0)), std(k1list(oocyteM == 1 & matEnv == 0)), .1);
set(ek10, 'color', 'k', 'linewidth', 2);

k20 = bar(.5, mean(k1list(oocyteM == 2 & matEnv == 0)), .25, 'facecolor', [.8 1 1]);
hold on;
ek20 = terrorbar(.5, mean(k1list(oocyteM == 2 & matEnv == 0)), std(k1list(oocyteM == 2 & matEnv == 0)), .1);
set(ek20, 'color', 'k', 'linewidth', 2);

% MM matured
k02 = bar(1, mean(k1list(oocyteM == 0 & matEnv == 2)), .25, 'facecolor', [.2 .6 .6]);
hold on;
ek02 = terrorbar(1, mean(k1list(oocyteM == 0 & matEnv == 2)), std(k1list(oocyteM == 0 & matEnv == 2)), .1);
set(ek02, 'color', 'k', 'linewidth', 2);

k12 = bar(1.25, mean(k1list(oocyteM == 1 & matEnv == 2)), .25, 'facecolor', [.4 .8 .8]);
hold on;
ek12 = terrorbar(1.25, mean(k1list(oocyteM == 1 & matEnv == 2)), std(k1list(oocyteM == 1 & matEnv == 2)), .1);
set(ek12, 'color', 'k', 'linewidth', 2);

k22 = bar(1.5, mean(k1list(oocyteM == 2 & matEnv == 2)), .25, 'facecolor', [.8 1 1]);
hold on;
ek22 = terrorbar(1.5, mean(k1list(oocyteM == 2 & matEnv == 2)), std(k1list(oocyteM == 2 & matEnv == 2)), .1);
set(ek22, 'color', 'k', 'linewidth', 2);

% KSOM matured
k01 = bar(2, mean(k1list(oocyteM == 0 & matEnv == 1)), .25, 'facecolor', [.2 .6 .6]);
hold on;
ek01 = terrorbar(2, mean(k1list(oocyteM == 0 & matEnv == 1)), std(k1list(oocyteM == 0 & matEnv == 1)), .1);
set(ek01, 'color', 'k', 'linewidth', 2);

k11 = bar(2.25, mean(k1list(oocyteM == 1 & matEnv == 1)), .25, 'facecolor', [.4 .8 .8]);
hold on;
ek11 = terrorbar(2.25, mean(k1list(oocyteM == 1 & matEnv == 1)), std(k1list(oocyteM == 1 & matEnv == 1)), .1);
set(ek11, 'color', 'k', 'linewidth', 2);

k21 = bar(2.5, mean(k1list(oocyteM == 2 & matEnv == 1)), .25, 'facecolor', [.8 1 1]);
hold on;
ek21 = terrorbar(2.5, mean(k1list(oocyteM == 2 & matEnv == 1)), std(k1list(oocyteM == 2 & matEnv == 1)), .1);
set(ek21, 'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [0.25 1.25 2.25])
set(gca, 'fontsize', 14);
ylim([0 .3]);
set(gca, 'xticklabel', {'in vivo', 'MM', 'KSOM'});
ylabel('k_1 parameter');
title('Effect of nuclear maturation on mechanics');
xlim([-0.5 3]);
%

legend([k01 k11 k21], {'GV', 'MI', 'MII'});

11111111
[p h] = ranksum(k1list(oocyteM == 2 & matEnv == 0), k1list(oocyteM == 2 & matEnv == 1))
[p h] = ranksum(k1list(oocyteM == 2 & matEnv == 0), k1list(oocyteM == 2 & matEnv == 2))
[p h] = ranksum(k1list(oocyteM == 2 & matEnv == 1), k1list(oocyteM == 2 & matEnv == 2))




%% Bar plot of MII oocytes matured in vivo, MM, or KSOM


figure(5);
clf;

% In Vivo matured
k20 = bar(.5, mean(k1list(oocyteM == 2 & matEnv == 0)), .8, 'facecolor', [.8 1 1]);
hold on;
ek20 = terrorbar(.5, mean(k1list(oocyteM == 2 & matEnv == 0)), std(k1list(oocyteM == 2 & matEnv == 0)), .1);
set(ek20, 'color', 'k', 'linewidth', 2);

% MM matured
k22 = bar(1.5, mean(k1list(oocyteM == 2 & matEnv == 2)), .8, 'facecolor', [.8 1 1]);
hold on;
ek22 = terrorbar(1.5, mean(k1list(oocyteM == 2 & matEnv == 2)), std(k1list(oocyteM == 2 & matEnv == 2)), .1);
set(ek22, 'color', 'k', 'linewidth', 2);

% KSOM matured
k21 = bar(2.5, mean(k1list(oocyteM == 2 & matEnv == 1)), .8, 'facecolor', [.8 1 1]);
hold on;
ek21 = terrorbar(2.5, mean(k1list(oocyteM == 2 & matEnv == 1)), std(k1list(oocyteM == 2 & matEnv == 1)), .1);
set(ek21, 'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [0.5 1.5 2.5])
set(gca, 'fontsize', 14);
ylim([0 .3]);
set(gca, 'xticklabel', {'in vivo', 'MM', 'KSOM'});
ylabel('k_1 parameter');
title('MII oocyte mechanics depend on environment');
xlim([0 3]);
%



%% Bar plot k1 of oocytes for non-fert/fert/blast


figure(5);
clf;
paramToPlot = n1list;


% Non-ferts
k00 = bar(.5, mean(paramToPlot(fertList == 0 & matEnv == 0)), .8, 'facecolor', [.85 .65 .2]);
hold on;
ek00 = terrorbar(.5, mean(paramToPlot(fertList == 0 & matEnv == 0)), std(paramToPlot(fertList == 0 & matEnv == 0)), .1);
set(ek00, 'color', 'k', 'linewidth', 2);

% Fert, no blast
k01 = bar(1.5, mean(paramToPlot(fertList == 1 & hatchList == 0 & matEnv == 0)), .8, 'facecolor', [.9 .2 .8]);
hold on;
ek01 = terrorbar(1.5, mean(paramToPlot(fertList == 1 & hatchList == 0 & matEnv == 0)), ...
    std(paramToPlot(fertList == 1 & hatchList == 0 & matEnv == 0)), .1);
set(ek01, 'color', 'k', 'linewidth', 2);

% Fert + hatched blast
k02 = bar(2.5, mean(paramToPlot(hatchList == 1 & matEnv == 0)), .8, 'facecolor', [0 .2 .8]);
hold on;
ek02 = terrorbar(2.5, mean(paramToPlot(hatchList == 1 & matEnv == 0)), ...
    std(paramToPlot(hatchList == 1 & matEnv == 0)), .1);
set(ek02, 'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [0.5 1.5 2.5])
set(gca, 'fontsize', 14);
ylim([0 12]);
set(gca, 'xticklabel', {'non-fert', 'fert', 'hatch'});
ylabel('k_1 parameter');
title('MII oocyte mechanics correlated with development');
xlim([0 3]);
%

[p h] = ranksum(paramToPlot(fertList == 0 & matEnv == 0), paramToPlot(fertList == 1 & hatchList == 0 & matEnv == 0))
[p h] = ranksum(paramToPlot(fertList == 0 & matEnv == 0), paramToPlot(hatchList == 1 & matEnv == 0))
[p h] = ranksum(paramToPlot(fertList == 1 & hatchList == 0 & matEnv == 0), paramToPlot(hatchList == 1 & matEnv == 0))





%% Plot proportion of MII over time


figure(7);
clf;

colorList = [[.2 .2 .6];[.8 .2 .8];[.6 .6 .6]];
p8 = sum(oocyteM == 2 & measHourList < 12 & matEnv == 0)/sum(measHourList < 12 & matEnv == 0);
p14 = sum(oocyteM == 2 & measHourList == 14 & matEnv == 0)/sum(measHourList == 14 & matEnv == 0);
p20 = sum(oocyteM == 2 & measHourList == 20 & matEnv == 0)/sum(measHourList == 20 & matEnv == 0);

se8 = sqrt(p8*(1-p8)/sum(measHourList < 12 & matEnv == 0));
se14 = sqrt(p14*(1-p14)/sum(measHourList == 14 & matEnv == 0));
se20 = sqrt(p20*(1-p20)/sum(measHourList == 20 & matEnv == 0));

% 12hr
k00 = bar(.5, p8, .8, 'facecolor', 1/4 + 3/4*colorList(1,:));
hold on;
ek00 = terrorbar(.5, p8, p8, 1.96*se8 + .5/sum(measHourList < 12), .1);
set(ek00, 'color', 'k', 'linewidth', 2);

% 16 hr
k01 = bar(1.5, p14, .8, 'facecolor', 1/4 + 3/4*colorList(2,:));
hold on;
ek01 = terrorbar(1.5, p14, 1.96*se14 + .5/sum(measHourList == 14), .1);
set(ek01, 'color', 'k', 'linewidth', 2);

% 20 hr
k02 = bar(2.5, p20, .8, 'facecolor', 1/4 + 3/4*colorList(3,:));
hold on;
ek02 = terrorbar(2.5, p20, 1.96*se20 + .5/sum(measHourList == 20), .1);
set(ek02, 'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [0.5 1.5 2.5])
set(gca, 'fontsize', 14);
ylim([0 1]);
set(gca, 'xticklabel', {'8 hr', '14 hr', '20 hr'});
ylabel('proportion at MII stage');
title('Morphology over time');
xlim([0 3]);
grid on;
%


%% Plot mechanical change over time (8hr, 14hr, 20hr)

figure(8);
clf;
hold on;
paramToPlot = n1list;

p8 = mean(paramToPlot(measHourList < 12 & matEnv == 0));
p14 = mean(paramToPlot(measHourList == 14 & matEnv == 0));
p20 = mean(paramToPlot(measHourList == 20 & matEnv == 0));

se8 = std(paramToPlot(measHourList < 12 & matEnv == 0));
se14 = std(paramToPlot(measHourList == 14 & matEnv == 0));
se20 = std(paramToPlot(measHourList == 20 & matEnv == 0));

sp8 = prctile(paramToPlot(measHourList < 12 & matEnv == 0), [10 25 50 75 90]);
sp14 = prctile(paramToPlot(measHourList == 14 & matEnv == 0), [10 25 50 75 90]);
sp20 = prctile(paramToPlot(measHourList == 20 & matEnv == 0), [10 25 50 75 90]);

xBar = NaN*zeros(1,length(paramToPlot));
xBar(measHourList < 12 & matEnv == 0) = 0;
xBar(measHourList == 14 & matEnv == 0) = 1;
xBar(measHourList == 20 & matEnv == 0) = 2;

% 8hr
% k00 = bar(.5, p8, .8, 'facecolor', 1/4 + 3/4*colorList(1,:));
%
% % 14 hr
% k01 = bar(1.5, p14, .8, 'facecolor', 1/4 + 3/4*colorList(2,:));
%
% % 20 hr
% k02 = bar(2.5, p20, .8, 'facecolor', 1/4 + 3/4*colorList(3,:));

% 8hr
plot([1 1], [sp8(1) sp8(5)], 'linewidth', 2, 'color', colorList(1,:)/2);
plot([.9 1.1], [sp8(5) sp8(5)], 'linewidth', 2, 'color', colorList(1,:)/2);
plot([.9 1.1], [sp8(1) sp8(1)], 'linewidth', 2, 'color', colorList(1,:)/2);
patch('xdata', [.7 .7 1.3 1.3], 'ydata', [sp8(2) sp8(4) sp8(4) sp8(2)], 'edgecolor', colorList(1,:)/2, ...
    'facecolor', 1/4 + 3/4*colorList(1,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
plot([.7 1.3], [sp8(3) sp8(3)], 'linewidth', 2, 'color', colorList(1,:)/2);

% 14hr
plot([2 2], [sp14(1) sp14(5)], 'linewidth', 2, 'color', colorList(2,:)/2);
plot([1.9 2.1], [sp14(5) sp14(5)], 'linewidth', 2, 'color', colorList(2,:)/2);
plot([1.9 2.1], [sp14(1) sp14(1)], 'linewidth', 2, 'color', colorList(2,:)/2);
patch('xdata', [1.7 1.7 2.3 2.3], 'ydata', [sp14(2) sp14(4) sp14(4) sp14(2)], 'edgecolor', colorList(2,:)/2, ...
    'facecolor', 1/4 + 3/4*colorList(2,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
plot([1.7 2.3], [sp14(3) sp14(3)], 'linewidth', 2, 'color', colorList(2,:)/2);

% 20hr
plot([3 3], [sp20(1) sp20(5)], 'linewidth', 2, 'color', colorList(3,:)/2);
plot([2.9 3.1], [sp20(5) sp20(5)], 'linewidth', 2, 'color', colorList(3,:)/2);
plot([2.9 3.1], [sp20(1) sp20(1)], 'linewidth', 2, 'color', colorList(3,:)/2);
patch('xdata', [2.7 2.7 3.3 3.3], 'ydata', [sp20(2) sp20(4) sp20(4) sp20(2)], 'edgecolor', colorList(3,:)/2, ...
    'facecolor', 1/4 + 3/4*colorList(3,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
plot([2.7 3.3], [sp20(3) sp20(3)], 'linewidth', 2, 'color', colorList(3,:)/2);


h = cell(1,length(k1list));
hLegend = cell(1,legendParamNums);
theta = 0:0.1:(2*pi);
xLocs = [1 2 3];
for i = 1:3
    
    currParam = [];
    
    if i == 1
        currParam = paramToPlot(measHourList < 12 & matEnv == 0);
    elseif i == 2
        currParam = paramToPlot(measHourList == 14 & matEnv == 0);
    else
        currParam = paramToPlot(measHourList == 20 & matEnv == 0);
    end
    
    for j = 1:length(currParam)
        xPatch = cos(theta)*.03 + xLocs(i) + rand(1)/2 - .25;
        yPatch = sin(theta)*.15 + currParam(j);
        
        h{i} = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', colorList(i,:)/2, ...
            'facecolor', colorList(i,:), 'facealpha', .3, 'edgealpha', .3);
    end
end


% kBoxPlot = boxplot(paramToPlot, xBar);

% ek00 = terrorbar(.5, p8, se8, .1);
% set(ek00, 'color', 'k', 'linewidth', 2);
% ek01 = terrorbar(1.5, p14, se14, .1);
% set(ek01, 'color', 'k', 'linewidth', 2);
% ek02 = terrorbar(2.5, p20, se20, .1);
% set(ek02, 'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [1 2 3])
set(gca, 'fontsize', 14);
% ylim([0 .2]);
ylim([0 11]);
set(gca, 'xticklabel', {'8 hr', '14 hr', '20 hr'});
ylabel('\eta_1 parameter');
title('Mechanics over time');
xlim([.5 3.5]);
grid on;


[p h] = ranksum(paramToPlot(measHourList < 12 & matEnv == 0), paramToPlot(measHourList == 14 & matEnv == 0))
[p h] = ranksum(paramToPlot(measHourList < 12 & matEnv == 0), paramToPlot(measHourList == 20 & matEnv == 0))
[p h] = ranksum(paramToPlot(measHourList == 14 & matEnv == 0), paramToPlot(measHourList == 20 & matEnv == 0))



%%  Plot mechanics and morphology over time (8, 14, 20 hrs)

figure(8);
clf;

% 8 hrs, GV
k00 = bar(0, mean(k1list(oocyteM == 0 & matEnv == 0 & measHourList < 12)), .25, 'facecolor', [.2 .6 .6]);
hold on;
ek00 = terrorbar(0, mean(k1list(oocyteM == 0 & matEnv == 0 & measHourList < 12)), ...
    std(k1list(oocyteM == 0 & matEnv == 0 & measHourList < 12)), .1);
set(ek00, 'color', 'k', 'linewidth', 2);

% 8 hrs, MI
k10 = bar(.25, mean(k1list(oocyteM == 1 & matEnv == 0 & measHourList < 12)), .25, 'facecolor', [.4 .8 .8]);
hold on;
ek10 = terrorbar(.25, mean(k1list(oocyteM == 1 & matEnv == 0 & measHourList < 12)), ...
    std(k1list(oocyteM == 1 & matEnv == 0 & measHourList < 12)), .1);
set(ek10, 'color', 'k', 'linewidth', 2);

% % 8 hrs, MII
% k20 = bar(.5, mean(k1list(oocyteM == 2 & matEnv == 0 & measHourList < 12)), .25, 'facecolor', [.8 1 1]);
% hold on;
% ek20 = terrorbar(.5, mean(k1list(oocyteM == 2 & matEnv == 0 & measHourList < 12)), ...
%     std(k1list(oocyteM == 2 & matEnv == 0 & measHourList < 12)), .1);
% set(ek20, 'color', 'k', 'linewidth', 2);

%%%%%%%%%%%%%%%%%

% 14 hrs, GV
k02 = bar(1, mean(k1list(oocyteM == 0 & matEnv == 0 & measHourList == 14)), .25, 'facecolor', [.2 .6 .6]);
hold on;
ek02 = terrorbar(1, mean(k1list(oocyteM == 0 & matEnv == 0 & measHourList == 14)), ...
    std(k1list(oocyteM == 0 & matEnv == 0 & measHourList == 14)), .1);
set(ek02, 'color', 'k', 'linewidth', 2);

% 14 hrs, MI
k12 = bar(1.25, mean(k1list(oocyteM == 1 & matEnv == 0 & measHourList == 14)), .25, 'facecolor', [.4 .8 .8]);
hold on;
ek12 = terrorbar(1.25, mean(k1list(oocyteM == 1 & matEnv == 0 & measHourList == 14)), ...
    std(k1list(oocyteM == 1 & matEnv == 0 & measHourList == 14)), .1);
set(ek12, 'color', 'k', 'linewidth', 2);

% 14 hrs, MII
k22 = bar(1.5, mean(k1list(oocyteM == 2 & matEnv == 0 & measHourList == 14)), .25, 'facecolor', [.8 1 1]);
hold on;
ek22 = terrorbar(1.5, mean(k1list(oocyteM == 2 & matEnv == 0 & measHourList == 14)), ...
    std(k1list(oocyteM == 2 & matEnv == 0 & measHourList == 14)), .1);
set(ek22, 'color', 'k', 'linewidth', 2);

%%%%%%%%%%%%%%%%%%

% 20 hrs, GV
% k01 = bar(2, mean(k1list(oocyteM == 0 & matEnv == 0 & measHourList == 20)), .25, 'facecolor', [.2 .6 .6]);
% hold on;
% ek01 = terrorbar(2, mean(k1list(oocyteM == 0 & matEnv == 0 & measHourList == 20)), ...
%     std(k1list(oocyteM == 0 & matEnv == 0 & measHourList == 20)), .1);
% set(ek01, 'color', 'k', 'linewidth', 2);

% 20 hrs, MI
k11 = bar(2.25, mean(k1list(oocyteM == 1 & matEnv == 0 & measHourList == 20)), .25, 'facecolor', [.4 .8 .8]);
hold on;
ek11 = terrorbar(2.25, mean(k1list(oocyteM == 1 & matEnv == 0 & measHourList == 20)), ...
    std(k1list(oocyteM == 1 & matEnv == 0 & measHourList == 20)), .1);
set(ek11, 'color', 'k', 'linewidth', 2);

% 20 hrs, MII
k21 = bar(2.5, mean(k1list(oocyteM == 2 & matEnv == 0 & measHourList == 20)), .25, 'facecolor', [.8 1 1]);
hold on;
ek21 = terrorbar(2.5, mean(k1list(oocyteM == 2 & matEnv == 0 & measHourList == 20)), ...
    std(k1list(oocyteM == 2 & matEnv == 0 & measHourList == 20)), .1);
set(ek21, 'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [0.25 1.25 2.25])
set(gca, 'fontsize', 14);
ylim([0 .25]);
set(gca, 'xticklabel', {'8 hrs', '14 hrs', '20 hrs'});
ylabel('k_1 parameter');
title('Effect of nuclear maturation on mechanics');
xlim([-0.5 3]);
%

legend([k00 k11 k21], {'GV', 'MI', 'MII'});

11111111

% are MIs and MIIs significant different at 14 hrs?
[p h] = ranksum(k1list(oocyteM == 2 & matEnv == 0 & measHourList == 14), ...
    k1list(oocyteM == 1 & matEnv == 0 & measHourList == 14))

% are MIIs at 14 hrs significant different from MIIs at 20 hrs?
[p h] = ranksum(k1list(oocyteM == 2 & matEnv == 0 & measHourList == 14), ...
    k1list(oocyteM == 2 & matEnv == 0 & measHourList == 20))



%% make scatterplot of fertilization vs mechanics

colorList = [[.86 .65 .2];[.2 .6 .8];[.8 .2 .8];[.2 .2 .6]];
h = cell(1,length(k1list));
theta = 0:0.1:(2*pi);

measHourFert = measHourList(~isnan(fertList) & matEnv == 0 & measHourList > 12);
kFert = k1list(~isnan(fertList) & matEnv == 0 & measHourList > 12);
nFert = n1list(~isnan(fertList) & matEnv == 0 &  measHourList > 12);
mFert = oocyteM(~isnan(fertList) & matEnv == 0 &  measHourList > 12);
outcomeFert1 = fertList(~isnan(fertList) & matEnv == 0 &  measHourList > 12);
outcomeFert2 = hatchList(~isnan(fertList) & matEnv == 0 &  measHourList > 12);

hLegend = cell(1,4);
legendParam = 0;


figure(1);
clf;
hold on;


for i = 1:length(kFert)
    
    if measHourFert(i) > 14 % late time point
        legendParam = 4;
        currColor = colorList(4,:);
    elseif outcomeFert2(i) == 1 % blast formed
        legendParam = 1;
        currColor = colorList(1,:);
    elseif outcomeFert1(i) == 1 % fertilized / PN formed
        legendParam = 2;
        currColor = colorList(2,:);
    else
        legendParam = 3;
        currColor = colorList(3,:);
    end
    
    
    xPatch = cos(theta)*.0013*2 + kFert(i);
    yPatch = sin(theta)*.17*2 + nFert(i);
    
    hLegend{legendParam} = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', currColor/2, ...
        'facecolor', currColor, 'facealpha', .5, 'edgealpha', .5);
    
    %         [xPatch, yPatch, zPatch] = sphere(10);
    %         xPatch = xPatch*.003 + kFert(i);
    %         yPatch = yPatch*.3 + nFert(i);
    %         zPatch = zPatch*.15 + mFert(i);
    %         currPoint = surf2patch(xPatch, yPatch, zPatch);
    %
    %         hLegend{legendParam} = patch(currPoint, 'edgecolor', currColor, ...
    %             'facecolor', currColor, 'facealpha', .3, 'edgealpha', 0);
end

% plot SVM decision boundary
axisLims = [.06 .145 -.5 9];
oocyteClassifier = fitcsvm([kFert; nFert]', outcomeFert1 > 0, 'ResponseName', ...
    'Viability', 'KernelFunction', 'rbf', 'KernelScale', .8, 'standardize', true);
[x1Grid,x2Grid] = meshgrid(linspace(axisLims(1),axisLims(2),100),...
    linspace(axisLims(3),axisLims(4),100));
xGrid = [x1Grid(:),x2Grid(:)];
[~,scores] = predict(oocyteClassifier,xGrid);
[c h] = contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),-.0*[1 1],'k', 'linewidth', 2);


set(gca, 'FontSize', 14);
title('Oocyte outcomes after IVF');
grid on;
axis([.06 .145 -.5 9]);
xlabel('k_1 parameter');
ylabel('\eta_1 parameter');
zlabel('morphology');

legend([hLegend{1}, hLegend{2}, hLegend{3} hLegend{4}], ...
    {'Blastocyst', 'Fertilized', 'No fertilization', '20hr time'}, 'Location', 'NorthEast');





%% Make subfigure for distance from decision boundary

[predictOut, scores] = predict(oocyteClassifier, [kFert; nFert]');
scores = scores(:,2);

% calculate min geometric distance 
xContourScaled = (c(1,2:end) - min(kFert))/(max(kFert) - min(kFert));
yContourScaled = (c(2,2:end) - min(nFert))/(max(nFert) - min(nFert));
xScaled = (kFert - min(kFert))/(max(kFert) - min(kFert));
yScaled = (nFert - min(nFert))/(max(nFert) - min(nFert));


decDistList = zeros(1,length(kFert));
for i = 1:length(decDistList)
    
    currDistList = zeros(1,length(xContourScaled));
    
    for j = 1:length(xContourScaled)
        currDistList(j) = sqrt((xScaled(i) - xContourScaled(j))^2 + ...
            (yScaled(i) - yContourScaled(j))^2);
    end
    
    decDistList(i) = min(currDistList);
end

decDistList = decDistList .* sign(scores)';


figure(10);
clf;
hold on;
paramToPlot = decDistList;

pBlast = paramToPlot(outcomeFert2 == 1 & measHourFert == 14);
pFert = paramToPlot(outcomeFert1 == 1 & outcomeFert2 == 0 & measHourFert == 14);
pNoFert = paramToPlot(outcomeFert1 == 0 & outcomeFert2 == 0 & measHourFert == 14);
pLate = paramToPlot(measHourFert > 14);

mBlast = mean(pBlast);
mFert = mean(pFert);
mNoFert = mean(pNoFert);
mLate = mean(pLate);

sBlast = std(pBlast);
sFert = std(pFert);
sNoFert = std(pNoFert);
sLate = std(pLate);

spBlast = prctile(pBlast, [10 25 50 75 90]);
spFert = prctile(pFert, [10 25 50 75 90]);
spNoFert = prctile(pNoFert, [10 25 50 75 90]);
spLate = prctile(pLate, [10 25 50 75 90]);

% Blast
% patch('xdata', [.7 .7 1.3 1.3 .7], 'ydata', [0 mBlast mBlast 0 0], ...
%     'linewidth', 2, 'edgecolor', colorList(1,:)/2, 'facecolor', 1/4 + 3/4*colorList(1,:), ...
%     'facealpha', 1, 'edgealpha', 1);
% plot([1 1], [mBlast-sBlast mBlast+sBlast], 'linewidth', 2, 'color', colorList(1,:)/2);
% plot([.9 1.1], [mBlast+sBlast mBlast+sBlast], 'linewidth', 2, 'color', colorList(1,:)/2);
% plot([.9 1.1], [mBlast-sBlast mBlast-sBlast], 'linewidth', 2, 'color', colorList(1,:)/2);

plot([1 1], [spBlast(1) spBlast(5)], 'linewidth', 2, 'color', colorList(1,:)/2);
plot([.9 1.1], [spBlast(5) spBlast(5)], 'linewidth', 2, 'color', colorList(1,:)/2);
plot([.9 1.1], [spBlast(1) spBlast(1)], 'linewidth', 2, 'color', colorList(1,:)/2);
patch('xdata', [.7 .7 1.3 1.3], 'ydata', [spBlast(2) spBlast(4) spBlast(4) spBlast(2)], 'edgecolor', colorList(1,:)/2, ...
    'facecolor', 1/3 + 2/3*colorList(1,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
plot([.7 1.3], [spBlast(3) spBlast(3)], 'linewidth', 2, 'color', colorList(1,:)/2);

% Fert
% patch('xdata', [1.7 1.7 2.3 2.3 1.7], 'ydata', [0 mFert mFert 0 0], ...
%     'linewidth', 2, 'edgecolor', colorList(2,:)/2, 'facecolor', 1/4 + 3/4*colorList(2,:), ...
%     'facealpha', 1, 'edgealpha', 1);
% plot([2 2], [mFert-sFert mFert+sFert], 'linewidth', 2, 'color', colorList(2,:)/2);
% plot([1.9 2.1], [mFert+sFert mFert+sFert], 'linewidth', 2, 'color', colorList(2,:)/2);
% plot([1.9 2.1], [mFert-sFert mFert-sFert], 'linewidth', 2, 'color', colorList(2,:)/2);

plot([2 2], [spFert(1) spFert(5)], 'linewidth', 2, 'color', colorList(2,:)/2);
plot([1.9 2.1], [spFert(5) spFert(5)], 'linewidth', 2, 'color', colorList(2,:)/2);
plot([1.9 2.1], [spFert(1) spFert(1)], 'linewidth', 2, 'color', colorList(2,:)/2);
patch('xdata', [1.7 1.7 2.3 2.3], 'ydata', [spFert(2) spFert(4) spFert(4) spFert(2)], 'edgecolor', colorList(2,:)/2, ...
    'facecolor', 1/3 + 2/3*colorList(2,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
plot([1.7 2.3], [spFert(3) spFert(3)], 'linewidth', 2, 'color', colorList(2,:)/2);

% No fert
% patch('xdata', [2.7 2.7 3.3 3.3 2.7], 'ydata', [0 mNoFert mNoFert 0 0], ...
%     'linewidth', 2, 'edgecolor', colorList(3,:)/2, 'facecolor', 1/4 + 3/4*colorList(3,:), ...
%     'facealpha', 1, 'edgealpha', 1);
% plot([3 3], [mNoFert-sNoFert mNoFert+sNoFert], 'linewidth', 2, 'color', colorList(3,:)/2);
% plot([2.9 3.1], [mNoFert+sNoFert mNoFert+sNoFert], 'linewidth', 2, 'color', colorList(3,:)/2);
% plot([2.9 3.1], [mNoFert-sNoFert mNoFert-sNoFert], 'linewidth', 2, 'color', colorList(3,:)/2);

plot([3 3], [spNoFert(1) spNoFert(5)], 'linewidth', 2, 'color', colorList(3,:)/2);
plot([2.9 3.1], [spNoFert(5) spNoFert(5)], 'linewidth', 2, 'color', colorList(3,:)/2);
plot([2.9 3.1], [spNoFert(1) spNoFert(1)], 'linewidth', 2, 'color', colorList(3,:)/2);
patch('xdata', [2.7 2.7 3.3 3.3], 'ydata', [spNoFert(2) spNoFert(4) spNoFert(4) spNoFert(2)], 'edgecolor', colorList(3,:)/2, ...
    'facecolor', 1/3 + 2/3*colorList(3,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
plot([2.7 3.3], [spNoFert(3) spNoFert(3)], 'linewidth', 2, 'color', colorList(3,:)/2);

% Late
% patch('xdata', [3.7 3.7 4.3 4.3 3.7], 'ydata', [0 mLate mLate 0 0], ...
%     'linewidth', 2, 'edgecolor', colorList(4,:)/2, 'facecolor', 1/4 + 3/4*colorList(4,:), ...
%     'facealpha', 1, 'edgealpha', 1);
% plot([4 4], [mLate-sLate mLate+sLate], 'linewidth', 2, 'color', colorList(4,:)/2);
% plot([3.9 4.1], [mLate+sLate mLate+sLate], 'linewidth', 2, 'color', colorList(4,:)/2);
% plot([3.9 4.1], [mLate-sLate mLate-sLate], 'linewidth', 2, 'color', colorList(4,:)/2);

plot([4 4], [spLate(1) spLate(5)], 'linewidth', 2, 'color', colorList(4,:)/2);
plot([3.9 4.1], [spLate(5) spLate(5)], 'linewidth', 2, 'color', colorList(4,:)/2);
plot([3.9 4.1], [spLate(1) spLate(1)], 'linewidth', 2, 'color', colorList(4,:)/2);
patch('xdata', [3.7 3.7 4.3 4.3], 'ydata', [spLate(2) spLate(4) spLate(4) spLate(2)], 'edgecolor', colorList(4,:)/2, ...
    'facecolor', 1/3 + 2/3*colorList(4,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
plot([3.7 4.3], [spLate(3) spLate(3)], 'linewidth', 2, 'color', colorList(4,:)/2);


hLegend = cell(1,4);
theta = 0:0.1:(2*pi);
xLocs = [1 2 3 4];

for i = 1:length(paramToPlot)
    
    if measHourFert(i) > 14 % late time point
        legendParam = 4;
        currColor = colorList(4,:);
    elseif outcomeFert2(i) == 1 % blast formed
        legendParam = 1;
        currColor = colorList(1,:);
    elseif outcomeFert1(i) == 1 && outcomeFert2(i) == 0 % fertilized / PN formed
        legendParam = 2;
        currColor = colorList(2,:);
    else
        legendParam = 3;
        currColor = colorList(3,:);
    end
    
    xPatch = cos(theta)*.04 + xLocs(legendParam) + rand(1)/2 - .25;
    yPatch = sin(theta)*.015 + paramToPlot(i);
    
    hLegend{legendParam} = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', currColor/2, ...
        'facecolor', currColor, 'facealpha', .4, 'edgealpha', .4);
    
end

line([0 5], [0 0], 'linewidth', 2', 'color', 'k', 'linestyle', '--');
set(gca, 'xtick', [1 2 3 4])
set(gca, 'fontsize', 14);
ylim([-.55 .8]);
set(gca, 'xticklabel', {'Blast', 'Fert', 'No fert', '20hr'});
ylabel('Distance from decision boundary');
title('Predicting fert and blast formation');
xlim([.5 4.5]);
grid on;



[p h] = ranksum(paramToPlot(outcomeFert2 == 1), paramToPlot(outcomeFert1 == 1 & outcomeFert2 == 0))
[p h] = ranksum(paramToPlot(outcomeFert1 == 1 & outcomeFert2 == 0), paramToPlot(outcomeFert1 == 0))
[p h] = ranksum(paramToPlot(outcomeFert1 == 0), paramToPlot(measHourFert > 14))



%% Plot proportion which fertilized vs controls


figure(11);
clf;

pMeasFert = sum(measHourList == 14 & fertList == 1)/sum(measHourList == 14 & ~isnan(fertList));
pMeasBlast = sum(measHourList == 14 & blastList == 1)/sum(measHourList == 14 & ~isnan(blastList) & fertList == 1);
pCtrlFert = 136/244;
pCtrlBlast = 86/136;

seMeasFert = sqrt(pMeasFert*(1-pMeasFert)/sum(measHourList == 14 & ~isnan(fertList)));
seMeasBlast = sqrt(pMeasBlast*(1-pMeasBlast)/sum(measHourList == 14 & ~isnan(blastList) & fertList == 1));
seCtrlFert = sqrt(pCtrlFert*(1-pCtrlFert)/244);
seCtrlBlast = sqrt(pCtrlBlast*(1-pCtrlBlast)/136);

% meas, fert
k00 = bar(.3, pMeasFert, .4, 'facecolor', 1/4 + 3/4*colorList(2,:));
hold on;
ek00 = terrorbar(.3, pMeasFert, 1.96*seMeasFert + .5/sum(measHourList == 14 & ~isnan(fertList)), .1);
set(ek00, 'color', 'k', 'linewidth', 2);

% ctrl, fert
k01 = bar(.7, pCtrlFert, .4, 'facecolor', [.6 .6 .6]);
hold on;
ek01 = terrorbar(.7, pCtrlFert, 1.96*seCtrlFert + .5/244, .1);
set(ek01, 'color', 'k', 'linewidth', 2);

% meas, blast
k10 = bar(1.3, pMeasBlast, .4, 'facecolor', 1/4 + 3/4*colorList(1,:));
hold on;
ek10 = terrorbar(1.3, pMeasBlast, 1.96*seMeasBlast + .5/sum(measHourList == 14 & ~isnan(blastList) & fertList == 1), .1);
set(ek10, 'color', 'k', 'linewidth', 2);

% ctrl, blast
k11 = bar(1.7, pCtrlBlast, .4, 'facecolor', [.6 .6 .6]);
hold on;
ek11 = terrorbar(1.7, pCtrlBlast, 1.96*seCtrlBlast + .5/136, .1);
set(ek11, 'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [.5 1.5])
set(gca, 'fontsize', 14);
set(gca, 'xticklabel', {'fertilized', 'blastocyst'});
ylim([0 1.2]);
ylabel('proportion');
title('fert and blast rates');
xlim([0 2]);
grid on;


%% Plot effect of maturation environment on mechanics

figure(2);
clf;
hold on;
paramToPlot = k1list;

colorList = [[.2 .2 .6];[.8 .2 .8];[.6 .6 .6]];
mCond = (oocyteM > 0);

pStart = paramToPlot(matEnv == 0 & measHourList < 12); % in vivo matured, 6 hr
pIVO_14 = paramToPlot(matEnv == 0 & measHourList == 14 & mCond);
pMM_14 = paramToPlot(matEnv == 2 & measHourList == 14 & mCond);
pKSOM_14 = paramToPlot(matEnv == 1 & measHourList == 14 & mCond);
pIVO_20 = paramToPlot(matEnv == 0 & measHourList == 20 & mCond);
pMM_20 = paramToPlot(matEnv == 2 & measHourList == 20 & mCond);
pKSOM_20 = paramToPlot(matEnv == 1 & measHourList == 20 & mCond);
pMM_26 = paramToPlot(matEnv == 2 & measHourList == 26 & mCond);
pKSOM_26 = paramToPlot(matEnv == 1 & measHourList == 26 & mCond);

spStart = prctile(pStart, [10 25 50 75 90]);
spIVO_14 = prctile(pIVO_14, [10 25 50 75 90]);
spMM_14 = prctile(pMM_14, [10 25 50 75 90]);
spKSOM_14 = prctile(pKSOM_14, [10 25 50 75 90]);
spIVO_20 = prctile(pIVO_20, [10 25 50 75 90]);
spMM_20 = prctile(pMM_20, [10 25 50 75 90]);
spKSOM_20 = prctile(pKSOM_20, [10 25 50 75 90]);
spMM_26 = prctile(pMM_26, [10 25 50 75 90]);
spKSOM_26 = prctile(pKSOM_26, [10 25 50 75 90]);

currParamList = {'spStart', 'spIVO_14', 'spMM_14', 'spKSOM_14', ...
    'spIVO_20', 'spMM_20', 'spKSOM_20', 'spMM_26', 'spKSOM_26'};
currXList = [1 3.2 4 4.8 7.2 8 8.8 11.1 11.9];
currNumList = [1 1 2 3 1 2 3 2 3];

hLegend = cell(1,3);
theta = 0:0.1:(2*pi);

% make box plot and scatter plot
for i = 1:length(currXList)
    
    currParam = eval(currParamList{i});
    fullParam = eval(currParamList{i}(2:end));
    currX = currXList(i);
    currNum = currNumList(i);
    plot([currX currX], [currParam(1) currParam(5)], 'linewidth', 2, 'color', colorList(currNum,:)/2);
    plot([currX-.1 currX+.1], [currParam(5) currParam(5)], 'linewidth', 2, 'color', colorList(currNum,:)/2);
    plot([currX-.1 currX+.1], [currParam(1) currParam(1)], 'linewidth', 2, 'color', colorList(currNum,:)/2);
    patch('xdata', [currX-.3 currX-.3 currX+.3 currX+.3], 'ydata', ...
        [currParam(2) currParam(4) currParam(4) currParam(2)], 'edgecolor', colorList(currNum,:)/2, ...
        'facecolor', 1/3 + 2/3*colorList(currNum,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
    plot([currX-.3 currX+.3], [currParam(3) currParam(3)], 'linewidth', 2, 'color', colorList(currNum,:)/2);

    
    for j = 1:length(fullParam)
        
        xPatch = cos(theta)*.05 + currX + rand(1)/2 - .25;
        yPatch = sin(theta)*.0012 + fullParam(j);
        
        hLegend{currNum} = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', colorList(currNum,:)/2, ...
            'facecolor', colorList(currNum,:), 'facealpha', .4, 'edgealpha', .4);
        
    end
    
end


line([0 12.5], [spStart(3) spStart(3)], 'linewidth', 2', 'color', 'k', 'linestyle', '--');

set(gca, 'xtick', [1 4 8 11.5])
set(gca, 'fontsize', 14);
ylim([0.05 .22]);
% ylim([-1 17]);
grid on;
set(gca, 'xticklabel', {'8hr', '14hr', '20hr', '26hr'});
legend([hLegend{1}, hLegend{2}, hLegend{3}], {'IVO', 'MM', 'KSOM'}, 'location', 'northwest');
ylabel('\eta_1 parameter');
title('Maturation environment affects oocyte mechanics');
xlim([0 12.5]);

[p h] = ranksum(eval(currParamList{1}(2:end)), eval(currParamList{5}(2:end)))
[p h] = ranksum(eval(currParamList{5}(2:end)), eval(currParamList{6}(2:end)))
[p h] = ranksum(eval(currParamList{5}(2:end)), eval(currParamList{7}(2:end)))
[p h] = ranksum(eval(currParamList{6}(2:end)), eval(currParamList{7}(2:end)))

[p h] = ranksum(eval(currParamList{8}(2:end)), eval(currParamList{9}(2:end)))
[p h] = ranksum(eval(currParamList{7}(2:end)), eval(currParamList{9}(2:end)))


%% Plot proportion of MII by environment over time


figure(3);
clf;
hold on;
paramToPlot = oocyteM;

colorList = [[.2 .2 .6];[.8 .2 .8];[.6 .6 .6]];

pStart = sum(matEnv == 0 & measHourList < 12 & paramToPlot == 2) / sum(matEnv == 0 & measHourList < 12); % in vivo matured, 6 hr
pIVO_14 = sum(matEnv == 0 & measHourList == 14 & paramToPlot == 2) / sum(matEnv == 0 & measHourList == 14);
pMM_14 = sum(matEnv == 2 & measHourList == 14 & paramToPlot == 2) / sum(matEnv == 2 & measHourList == 14);
pKSOM_14 = sum(matEnv == 1 & measHourList == 14 & paramToPlot == 2) / sum(matEnv == 1 & measHourList == 14);
pIVO_20 = sum(matEnv == 0 & measHourList == 20 & paramToPlot == 2) / sum(matEnv == 0 & measHourList == 20);
pMM_20 = sum(matEnv == 2 & measHourList == 20 & paramToPlot == 2) / sum(matEnv == 2 & measHourList  == 20);
pKSOM_20 = sum(matEnv == 1 & measHourList == 20 & paramToPlot == 2) / sum(matEnv == 1 & measHourList  == 20);
pMM_26 = sum(matEnv == 2 & measHourList == 26 & paramToPlot == 2) / sum(matEnv == 2 & measHourList  == 26);
pKSOM_26 = sum(matEnv == 1 & measHourList == 26 & paramToPlot == 2) / sum(matEnv == 1 & measHourList  == 26);

spStart = sqrt(pStart*(1-pStart)/sum(matEnv == 0 & measHourList < 12));
spIVO_14 = sqrt(pIVO_14*(1-pIVO_14)/sum(matEnv == 0 & measHourList == 14));
spMM_14 = sqrt(pMM_14*(1-pMM_14)/sum(matEnv == 2 & measHourList == 14));
spKSOM_14 = sqrt(pKSOM_14*(1-pKSOM_14)/sum(matEnv == 1 & measHourList == 14));
spIVO_20 = sqrt(pIVO_20*(1-pIVO_20)/sum(matEnv == 0 & measHourList == 20));
spMM_20 = sqrt(pMM_20*(1-pMM_20)/sum(matEnv == 2 & measHourList == 20));
spKSOM_20 = sqrt(pKSOM_20*(1-pKSOM_20)/sum(matEnv == 1 & measHourList == 20));
spMM_26 = sqrt(pMM_26*(1-pMM_26)/sum(matEnv == 2 & measHourList == 26));
spKSOM_26 = sqrt(pKSOM_26*(1-pKSOM_26)/sum(matEnv == 1 & measHourList == 26));


currParamList = {'pStart', 'pIVO_14', 'pMM_14', 'pKSOM_14', ...
    'pIVO_20', 'pMM_20', 'pKSOM_20', 'pMM_26', 'pKSOM_26'};
currXList = [1 3.2 4 4.8 7.2 8 8.8 11.1 11.9];
currNumList = [1 1 2 3 1 2 3 2 3];


hLegend = cell(1,3);
theta = 0:0.1:(2*pi);

% make box plot and scatter plot
for i = 1:length(currXList)
    
    currParam = eval(currParamList{i});
    currParamS = eval(['s' currParamList{i}]);
    currX = currXList(i);
    currNum = currNumList(i);

    
    
    hLegend{currNum} = patch('xdata', [currX-.3 currX-.3 currX+.3 currX+.3], 'ydata', ...
        [0 currParam currParam 0], 'edgecolor', colorList(currNum,:)/2, ...
        'facecolor', 1/3 + 2/3*colorList(currNum,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
    
    plot([currX-.1 currX+.1], (currParam + currParamS)*[1 1], ...
        'linewidth', 2, 'color', colorList(currNum,:)/2);
    plot([currX-.1 currX+.1], (currParam - currParamS)*[1 1], ...
        'linewidth', 2, 'color', colorList(currNum,:)/2);
    plot([currX currX], currParam + currParamS*[-1 1], 'linewidth', 2, 'color', colorList(currNum,:)/2);
    
end

set(gca, 'xtick', [1 4 8 11.5])
set(gca, 'fontsize', 14);
ylim([0 1.2]);
grid on;
set(gca, 'xticklabel', {'8hr', '14hr', '20hr', '26hr'});
legend([hLegend{1}, hLegend{2}, hLegend{3}], {'IVO', 'MM', 'KSOM'}, 'location', 'northeast');
ylabel('proportion at MII stage');
title('IVM morphology');
xlim([0 12.5]);


%% Break down by GV, MI, MII (early), MII (late)

% cat 1: GV 8-14 hr
% cat 2: MI 8-14 hr
% cat 3: MII 8-14 hr
% cat 4: MII 20 hr
% all in vivo matured


figure(4);
clf;
hold on;
paramToPlot = k1list;

colorList = [[.2 .2 .6];[.8 .2 .8];[.6 .6 .6]; [.3 .3 .3]];

pGV = paramToPlot(matEnv == 0 & measHourList < 20 & ~isnan(paramToPlot) & oocyteM == 0);
pMI = paramToPlot(matEnv == 0 & measHourList < 20 & ~isnan(paramToPlot) & oocyteM == 1);
pMIIearly = paramToPlot(matEnv == 0 & measHourList < 20 & ~isnan(paramToPlot) & oocyteM == 2);
pMIIlate = paramToPlot(matEnv == 0 & measHourList == 20 & ~isnan(paramToPlot) & oocyteM == 2);

meanList = [mean(pGV) mean(pMI) mean(pMIIearly) mean(pMIIlate)];
stdList = [std(pGV) std(pMI) std(pMIIearly) std(pMIIlate)];
currXList = [1 2 3 5];
currNumList = [1 2 3 4];

hLegend = cell(1,4);
theta = 0:0.1:(2*pi);

% make box plot and scatter plot
for i = 1:length(currXList)
    
    currParam = meanList(i);
    currParamS = stdList(i);
    currX = currXList(i);
    currNum = currNumList(i);
    
    hLegend{currNum} = patch('xdata', [currX-.3 currX-.3 currX+.3 currX+.3], 'ydata', ...
        [0 currParam currParam 0], 'edgecolor', colorList(currNum,:)/2, ...
        'facecolor', 1/3 + 2/3*colorList(currNum,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
    
    plot([currX-.1 currX+.1], (currParam + currParamS)*[1 1], ...
        'linewidth', 2, 'color', colorList(currNum,:)/2);
    plot([currX-.1 currX+.1], (currParam - currParamS)*[1 1], ...
        'linewidth', 2, 'color', colorList(currNum,:)/2);
    plot([currX currX], currParam + currParamS*[-1 1], 'linewidth', 2, 'color', colorList(currNum,:)/2);
    
end

set(gca, 'xtick', [1 2 3 5])
set(gca, 'fontsize', 14);
ylim([0 .2]);
grid on;
set(gca, 'xticklabel', {'GV', 'MI', 'MII', 'MII'});
% legend([hLegend{1}, hLegend{2}, hLegend{3}], {'IVO', 'MM', 'KSOM'}, 'location', 'northeast');
ylabel('k_1 parameter');
title('Mechanics vs morphology');
xlim([0 6]);

[h p] = ranksum(pGV, pMI)
[h p] = ranksum(pMIIearly, pMI)


%% Plot k1 for each category


% cat 1: GV 8-14 hr
% cat 2: MI 8-14 hr
% cat 3: MII 8-14 hr
% cat 4: MII 20 hr
% all in vivo matured


figure(4);
clf;
hold on;
paramToPlot = k1list;

colorList = [[.6 .6 .6];[.86 .65 .2];[.2 .6 .8];[.8 .2 .8];[.2 .2 .6]];

p8 = paramToPlot(matEnv == 0 & measHourList < 14 & ~isnan(paramToPlot));
p14B = paramToPlot(matEnv == 0 & measHourList == 14 & blastList == 1);
p14F = paramToPlot(matEnv == 0 & measHourList == 14 & fertList == 1);% & blastList == 0);
p14N = paramToPlot(matEnv == 0 & measHourList == 14 & fertList == 0);
p20 = paramToPlot(matEnv == 0 & measHourList == 20 & ~isnan(paramToPlot));

meanList = [mean(p8) mean(p14B) mean(p14F) mean(p14N) mean(p20)];
stdList = [std(p8) std(p14B) std(p14F) std(p14N) std(p20)];
currXList = [1 3 4 5 7];
currNumList = [1 2 3 4 1];

hLegend = cell(1,5);
theta = 0:0.1:(2*pi);

% make box plot and scatter plot
for i = 1:length(currXList)
    
    currParam = meanList(i);
    currParamS = stdList(i);
    currX = currXList(i);
    currNum = currNumList(i);
    
    hLegend{currNum} = patch('xdata', [currX-.3 currX-.3 currX+.3 currX+.3], 'ydata', ...
        [0 currParam currParam 0], 'edgecolor', colorList(currNum,:)/2, ...
        'facecolor', 1/3 + 2/3*colorList(currNum,:), 'facealpha', 1, 'edgealpha', 1, 'linewidth', 2);
    
    plot([currX-.1 currX+.1], (currParam + currParamS)*[1 1], ...
        'linewidth', 2, 'color', colorList(currNum,:)/2);
    plot([currX-.1 currX+.1], (currParam - currParamS)*[1 1], ...
        'linewidth', 2, 'color', colorList(currNum,:)/2);
    plot([currX currX], currParam + currParamS*[-1 1], 'linewidth', 2, 'color', colorList(currNum,:)/2);
    
end

set(gca, 'xtick', [1 4 7])
set(gca, 'fontsize', 14);
ylim([0 .22]);
grid on;
set(gca, 'xticklabel', {'8hr', '14hr', '20hr'});
legend([hLegend{2}, hLegend{3}, hLegend{4}], ...
    {'blast', 'fert','no fert'}, 'location', 'northwest');
ylabel('k_1 parameter');
title('Mechanics vs morphology');
xlim([0 8]);

[h p] = ranksum(p14B, p14F)
[h p] = ranksum(p14F, p14N)
[h p] = ranksum(p14B, p14N)

p14FNB = paramToPlot(matEnv == 0 & measHourList == 14 & fertList == 1 & blastList == 0);
p14NB = paramToPlot(matEnv == 0 & measHourList == 14 & blastList == 0);


[h p] = ranksum(p14F, p14N)
[h p] = ranksum(p14B, p14NB)

















