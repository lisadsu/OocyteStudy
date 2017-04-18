% Livia Zarnescu
% 10-7-2015
% Make some comparisons between zona, membrane, and bulk props

clear all;
close all;

%% 1. Compare bulk properties at different pressures

% 10-2-15 data
date = '10-2-15';
eNums = [4 5];
pressureVals = [1 2 3 4];
baseDir = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte';
dataDir = [baseDir '\' date ' analysis\AutoMeasure'];
colorVec = [0 0 .6; 0 .6 0; .6 0 0; 0 .6 .6];
markerList = {'+', 'o', '^', '*'};

k1M = zeros(length(eNums),length(pressureVals));
n1M = zeros(length(eNums),length(pressureVals));
k0M = zeros(length(eNums),length(pressureVals));
tauM = zeros(length(eNums),length(pressureVals));

f1 = figure(1);
clf;

% plot aspiration curves for different pressures
for i = 1:length(eNums)
    for j = 1:length(pressureVals)
        
        currFile = [dataDir, '\aspiration_data_', strrep(date, '-', '_'), ...
            '_E', num2str(eNums(i)), '_pt', num2str(pressureVals(j)), '.mat'];
        
        load(currFile);
        %%%%%%%%%%%%%%%%%%%
                
        pressureApplied = pressureVals(j)*0.1*0.85;
        Fin = pressureApplied * 6895 * pi * (20*128/126*10^-6)^2;
        [xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 1, [k0 k1 tau n1]');
        pause;
        %%%%%%%%%%%%%%%%%%%%

        k1M(i,j) = k1;
        n1M(i,j) = n1;
        k0M(i,j) = k0;
        tauM(i,j) = tau;
        
        figure(1);
        hold on;
        scatter(t, A*10^6, 150, 'marker', markerList{j}, 'linewidth', 2, ...
            'markeredgecolor', colorVec(i,:));
        plot(xfine-min(xfine), yfit*10^6, 'linewidth', 2, 'color', colorVec(i,:));
        
    end
end

% plot how parameters change for different pressures
f2 = figure(2);
clf;

subplot(2,2,1);
set(gca, 'fontsize', 14);
bar([.1 .2 .3 .4], k1M')
colormap([0 0 .6; 0 .6 0]);
title('k1');
xlim([0 .5]);

subplot(2,2,2);
set(gca, 'fontsize', 14);
bar([.1 .2 .3 .4], n1M')
colormap([0 0 .6; 0 .6 0]);
title('n1');
xlim([0 .5]);

subplot(2,2,3);
set(gca, 'fontsize', 14);
bar([.1 .2 .3 .4], k0M')
colormap([0 0 .6; 0 .6 0]);
title('k0');
xlim([0 .5]);

subplot(2,2,4);
set(gca, 'fontsize', 14);
bar([.1 .2 .3 .4], tauM')
colormap([0 0 .6; 0 .6 0]);
title('tau');
xlim([0 .5]);


%% 2. correlate zona stiffness and bulk props

dateList = {'10-8-15', '10-21-15'};
numEmbryos = 0;
mList = [];
eNums = [];
eFileNames = [];
dates = [];
extraString = [];
baseDir = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte';
colorVec = [0 .6 0; 0 0 .6; .6 0 0; 0 .6 .6];
markerList = {'+', 'o', '^', '*'};

% 10-8-15 data
% E3, E12 fragmented
% mList is for E1-20 at 24 hrs
mList = [mList 1 3 NaN 1 1 1 1 1 1 2 2 NaN 1 1 2 2 1 2 2 2]; % 1 = M2, 2 = M1, 3 = GV
numEmbryos = numEmbryos + 20;
eNums = [eNums 1:20];
eFileNames = [eFileNames 1:20];
dates = [dates ones(1,20)];
extraString = [extraString ones(1,20)];

% 10-21-15
mList = [mList 3 3 2 3 3]; % 1 = M2, 2 = M1, 3 = GV
numEmbryos = numEmbryos + 5;
eNums = [eNums 1:5];
eFileNames = [eFileNames 131:135];
dates = [dates 2*ones(1,5)];
extraString = [extraString zeros(1,5)];

colorList = zeros(numEmbryos, 3);
k1list = zeros(1,numEmbryos);
n1list = zeros(1,numEmbryos);
k0list = zeros(1,numEmbryos);
taulist = zeros(1,numEmbryos);

% tensionList = zeros(1,length(eNums));
zonaE = zeros(1,numEmbryos);
zonaDepth = zeros(1,numEmbryos); % depth @ .2 psi

% build vectors for each embryo
for i = 1:numEmbryos
    
    if ~isnan(mList(i))
        colorList(i,:) = colorVec(mList(i),:); %[0 0 .5]; %
        %         if (i < 6) || (i > 14)
        %             colorList(i,:) = colorList(i,:) + .3;
        %         end
    end
    
    extraStringCurr = '';
    if extraString(i)
        extraStringCurr = '_24';
    end
    
    % bulk mech data first
    currDate = dateList{dates(i)};
    dataDir = [baseDir '\' currDate ' analysis'];
    dataBulk = [dataDir, '\AutoMeasure\aspiration_data_', ...
        strrep(currDate, '-', '_'), ...
        '_E', num2str(eFileNames(i)), extraStringCurr, '.mat'];
    
    if ~exist(dataBulk, 'file')
        k1list(i) = NaN;
        n1list(i) = NaN;
        k0list(i) = NaN;
        taulist(i) = NaN;
    else
        load(dataBulk);
        k1list(i) = k1;
        n1list(i) = n1;
        k0list(i) = k0;
        taulist(i) = tau;
    end
    
    % surf tension currently all in one file, maybe add later
    
    % zona only data
    dataZona = [dataDir, '\Zona\allZona.mat'];
    
    if ~exist(dataZona, 'file')
        zonaE(i) = NaN;
        zonaDepth(i) = NaN;
    else
        
        load(dataZona);
        
        % if the cell was measured
        if ~isempty(pressureAll{eNums(i)})
            pSub = pressureAll{eNums(i)} - min(pressureAll{eNums(i)});
            dSub = depthAll{eNums(i)} - min(depthAll{eNums(i)});
            
            
            
            dInterp = 0:.1:10; % only look at small deformations (less than zona thickness)
            pInterp = interp1(dSub, pSub, dInterp, 'spline', 'extrap');
            p = polyfit(dInterp, pInterp, 1);
            
            % E = dP/dL * L0 = p(1) * L0
            zonaE(i) = cellSize(eNums(i))*(p(1)*6895); % units Pascals
            zonaDepth(i) = interp1(pSub, dSub, .2, 'spline', 'extrap'); % units um
            
        else
            zonaE(i) = NaN;
            zonaDepth(i) = NaN;
        end
        
    end
    
    
end


% now plot
figure(3);
clf;

% zonaE = zonaDepth*1000
cellsToPlot = ~isnan(zonaE) & ~isnan(k1list);

subplot(2,2,1);
set(gca, 'fontsize', 14);
scatter(k1list, zonaE/1000, 150, colorList, 'linewidth', 2, 'marker', '+');
ylabel('Zona Young''s modulus (kPa)');
xlabel('k1 bulk');
[r p] = corrcoef(k1list(cellsToPlot), zonaE(cellsToPlot)/1000)
title(['r^2 = ' num2str(r(1,2)^2)]);
grid on;

subplot(2,2,2);
set(gca, 'fontsize', 14);
scatter(n1list, zonaE/1000, 150, colorList, 'linewidth', 2, 'marker', '+');
ylabel('Zona Young''s modulus (kPa)');
xlabel('n1 bulk');
[r p] = corrcoef(n1list(cellsToPlot), zonaE(cellsToPlot)/1000)
title(['r^2 = ' num2str(r(1,2)^2)]);
grid on;

subplot(2,2,3);
set(gca, 'fontsize', 14);
scatter(k0list, zonaE/1000, 150, colorList, 'linewidth', 2, 'marker', '+');
ylabel('Zona Young''s modulus (kPa)');
xlabel('k0 bulk');
[r p] = corrcoef(k0list(cellsToPlot), zonaE(cellsToPlot)/1000)
title(['r^2 = ' num2str(r(1,2)^2)]);
grid on;

subplot(2,2,4);
set(gca, 'fontsize', 14);
scatter(taulist, zonaE/1000, 150, colorList, 'linewidth', 2, 'marker', '+');
ylabel('Zona Young''s modulus (kPa)');
xlabel('tau bulk');
[r p] = corrcoef(taulist(cellsToPlot), zonaE(cellsToPlot)/1000)
title(['r^2 = ' num2str(r(1,2)^2)]);
grid on;


%% 2.5. Correlate surface tension with bulk props

dateList = {'10-21-15'};
numEmbryos = 0;
mList = [];
eNums = [];
eFileNames = [];
dates = [];
extraString = [];
baseDir = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte';
colorVec = [.6 0 0; 0 0 .6; 0 .6 0; 0 .6 .6];
markerList = {'+', 'o', '^', '*'};

% 10-21-15 KSOM
mList = [mList 1 3 NaN 2 2 1 3 2 2 3 ]; % 1 = GV, 2 = M1, 3 = M2
numEmbryos = numEmbryos + 10;
eNums = [eNums 1:10];
eFileNames = [eFileNames 1:10];
dates = [dates ones(1,10)];
extraString = [extraString ones(1,10)];

% 10-21-15 MM
mList = [mList 3 1 NaN 3 3 3 3 3 3 3];
numEmbryos = numEmbryos + 10;
eNums = [eNums 11:20];
eFileNames = [eFileNames 21:30];
dates = [dates ones(1,10)];
extraString = [extraString ones(1,10)];

colorList = zeros(numEmbryos, 3);
k1list = zeros(1,numEmbryos);
n1list = zeros(1,numEmbryos);
k0list = zeros(1,numEmbryos);
taulist = zeros(1,numEmbryos);

% tensionList = zeros(1,length(eNums));
surfTensionAll = zeros(1,numEmbryos);
eqPressureAll = zeros(1,numEmbryos); % pressure for depth = 20 um

% build vectors for each embryo
for i = 1:numEmbryos
    
    if ~isnan(mList(i))
        colorList(i,:) = colorVec(mList(i),:); %[0 0 .5]; %
    end
    
    extraStringCurr = '';
    if extraString(i)
        extraStringCurr = '_24';
    end
    
    % bulk mech data first
    currDate = dateList{dates(i)};
    dataDir = [baseDir '\' currDate ' analysis'];
    dataBulk = [dataDir, '\AutoMeasure\aspiration_data_', ...
        strrep(currDate, '-', '_'), ...
        '_E', num2str(eFileNames(i)), extraStringCurr, '.mat'];
    
    if ~exist(dataBulk, 'file')
        k1list(i) = NaN;
        n1list(i) = NaN;
        k0list(i) = NaN;
        taulist(i) = NaN;
    else
        load(dataBulk);
        k1list(i) = k1;
        n1list(i) = n1;
        k0list(i) = k0;
        taulist(i) = tau;
    end
    

    % surf tension currently all in one file, maybe add later
    dataMembrane = [dataDir, '\Membrane\allMembrane.mat'];
    
    if ~exist(dataMembrane, 'file')
        surfTensionAll(i) = NaN;
        eqPressureAll(i) = NaN;
    else
        
        load(dataMembrane);
        
        if length(surfTension) >= eNums(i) && surfTension(eNums(i)) > 0
            surfTensionAll(i) = surfTension(eNums(i));
        else
            surfTensionAll(i) = NaN;
        end
        
        if length(eqPressure) >= eNums(i) && eqPressure(eNums(i)) > 0
            eqPressureAll(i) = eqPressure(eNums(i));
        else
            eqPressureAll(i) = NaN;
        end
        
    end
end


% now plot
figure(7);
clf;

% zonaE = zonaDepth*1000
cellsToPlot = ~isnan(surfTensionAll) & ~isnan(k1list);
paramToPlot = surfTensionAll;

n1list = log(n1list);
k0list = log(k0list);

subplot(2,2,1);
set(gca, 'fontsize', 14);
scatter(k1list(cellsToPlot), paramToPlot(cellsToPlot), ...
    150, colorList(cellsToPlot,:), 'linewidth', 2, 'marker', '+');
ylabel('Cell surface tension (N/m)');
xlabel('k1 bulk');
[r p] = corrcoef(k1list(cellsToPlot), paramToPlot(cellsToPlot))
title(['r^2 = ' num2str(r(1,2)^2)]);
grid on;

subplot(2,2,2);
set(gca, 'fontsize', 14);
scatter(n1list(cellsToPlot), paramToPlot(cellsToPlot), ...
    150, colorList(cellsToPlot,:), 'linewidth', 2, 'marker', '+');
ylabel('Cell surface tension (N/m)');
xlabel('n1 bulk');
[r p] = corrcoef(n1list(cellsToPlot), paramToPlot(cellsToPlot))
title(['r^2 = ' num2str(r(1,2)^2)]);
grid on;

subplot(2,2,3);
set(gca, 'fontsize', 14);
scatter(k0list(cellsToPlot), paramToPlot(cellsToPlot), ...
    150, colorList(cellsToPlot,:), 'linewidth', 2, 'marker', '+');
ylabel('Cell surface tension (N/m)');
xlabel('k0 bulk');
[r p] = corrcoef(k0list(cellsToPlot), paramToPlot(cellsToPlot))
title(['r^2 = ' num2str(r(1,2)^2)]);
grid on;

subplot(2,2,4);
set(gca, 'fontsize', 14);
scatter(taulist(cellsToPlot), paramToPlot(cellsToPlot), ...
    150, colorList(cellsToPlot,:), 'linewidth', 2, 'marker', '+');
ylabel('Cell surface tension (N/m)');
xlabel('tau bulk');
[r p] = corrcoef(taulist(cellsToPlot), paramToPlot(cellsToPlot))
title(['r^2 = ' num2str(r(1,2)^2)]);
grid on;


%% 3. compare different culture conditions (KSOM vs culture media)
% as well as time inside/outside incubator and mechanical measurement
% mech group was measured at
% exp group was measured only once at 24 hrs, but brought outside incubator
%   while mech group was being measured
% inc group was kept inside incubator and only measured once at 24 hrs

baseDir = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte';
groups = {'mechK', 'mechM', 'expK', 'expM', 'incK', 'incM'};
dates = {'10-8-15', '10-21-15'};
groupAll = [];
dateAll = [];
embryoNumAll = [];

% 10-8-15 data
% which index of groups is embryo in?
groupAll = [groupAll 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 5*ones(1,16) ...
    6*ones(1,12) 4*ones(1,9) 3*ones(1,8)];
dateAll = [dateAll ones(1,65)];
embryoNumAll = [embryoNumAll 1:65];

% 10-21-15
groupAll = [groupAll ones(1,10) 2*ones(1,20) 3*ones(1,11) ...
    4*ones(1,10) 6*ones(1,10) 5*ones(1,10)];
dateAll = [dateAll 2*ones(1,71)];
embryoNumAll = [embryoNumAll 1:71];

% mech is blue, exp is orange, inc is purple
colorVec = [0 0 .5; .6 .6 .8; 0 .5 .5; .6 .8 .8; 0 .5 0; .6 .8 .6];

numEmbryos = length(groupAll);
k1list = zeros(1,numEmbryos);
n1list = zeros(1,numEmbryos);
k0list = zeros(1,numEmbryos);
taulist = zeros(1,numEmbryos);

% build vectors for each embryo
for i = 1:numEmbryos
    
    if ~isnan(groupAll(i))
        colorList(i,:) = colorVec(groupAll(i),:);
    end
    
    % bulk mech data first
    extraString = '';
    if groupAll(i) < 3
        extraString = '_24';
    end
    
    dataDir = [baseDir '\' dates{dateAll(i)} ' analysis'];
    dataBulk = [dataDir, '\AutoMeasure\aspiration_data_', ...
        strrep(dates{dateAll(i)}, '-', '_'), ...
        '_E', num2str(embryoNumAll(i)), extraString, '.mat'];
    
    if ~exist(dataBulk, 'file')
        k1list(i) = NaN;
        n1list(i) = NaN;
        k0list(i) = NaN;
        taulist(i) = NaN;
    else
        load(dataBulk);
        k1list(i) = k1;
        n1list(i) = n1;
        k0list(i) = k0;
        taulist(i) = tau;
    end
end

% now bar graph comparing mech, exp, and inc only
k1Mech = k1list((groupAll < 3) & ~isnan(k1list));
k1Exp = k1list((groupAll > 2) & (groupAll < 5) & ~isnan(k1list));
k1Inc = k1list((groupAll > 4) & ~isnan(k1list));
k1K = k1list((mod(groupAll, 2) == 1) & ~isnan(k1list));
k1M = k1list((mod(groupAll, 2) == 0) & ~isnan(k1list));

figure(1);
clf;
set(gca, 'fontsize', 14);
h1 = bar(1, mean(k1Mech), .8, 'facecolor', colorVec(1,:));
hold on;
e1 = errorbar(1, mean(k1Mech), std(k1Mech),'color', 'k', 'linewidth', 2);

h2 = bar(2, mean(k1Exp), .8, 'facecolor', colorVec(3,:));
hold on;
e2 = errorbar(2, mean(k1Exp), std(k1Exp),'color', 'k', 'linewidth', 2);

h3 = bar(3, mean(k1Inc), .8, 'facecolor', colorVec(5,:));
hold on;
e3 = errorbar(3, mean(k1Inc), std(k1Inc),'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [1 2 3])
legend([h1 h2 h3], {'mech + outside incubator', 'outside incubator', ...
    'inside incubator'}, 'location', 'northwest');
set(gca, 'xticklabel', {'', '', ''});
xlim([.5 3.5]);
ylim([0 .22]);
grid on;
ylabel('stiffness ');
title('Oocyte stiffness after mechanical measurement');


figure(2);
clf;
set(gca, 'fontsize', 14);
h1 = bar(1, mean(k1K), .8, 'facecolor', colorVec(1,:));
hold on;
e1 = errorbar(1, mean(k1K), std(k1K),'color', 'k', 'linewidth', 2);

h2 = bar(2, mean(k1M), .8, 'facecolor', colorVec(2,:));
hold on;
e2 = errorbar(2, mean(k1M), std(k1M),'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [1 2])
legend([h1 h2], {'KSOM', 'MM'}, 'location', 'northwest');
set(gca, 'xticklabel', {'', ''});
xlim([.5 2.5]);
ylim([0 .2]);
grid on;
ylabel('stiffness ');
title('Oocyte stiffness depends on maturation conditions');

% allSeparated!
k11 = k1list((groupAll == 1) & ~isnan(k1list));
k12 = k1list((groupAll == 2) & ~isnan(k1list));
k13 = k1list((groupAll == 3) & ~isnan(k1list));
k14 = k1list((groupAll == 4) & ~isnan(k1list));
k15 = k1list((groupAll == 5) & ~isnan(k1list));
k16 = k1list((groupAll == 6) & ~isnan(k1list));

figure(3);
clf;
set(gca, 'fontsize', 14);
h1 = bar(.8, mean(k11), .4, 'facecolor', colorVec(1,:));
hold on;
e1 = errorbar(.8, mean(k11), std(k11), 'color', 'k', 'linewidth', 2);

h2 = bar(1.2, mean(k12), .4, 'facecolor', colorVec(2,:));
hold on;
e2 = errorbar(1.2, mean(k12), std(k12), 'color', 'k', 'linewidth', 2);

h3 = bar(1.8, mean(k13), .4, 'facecolor', colorVec(3,:));
hold on;
e3 = errorbar(1.8, mean(k13), std(k13), 'color', 'k', 'linewidth', 2);

h4 = bar(2.2, mean(k14), .4, 'facecolor', colorVec(4,:));
hold on;
e3 = errorbar(2.2, mean(k14), std(k14), 'color', 'k', 'linewidth', 2);

h5 = bar(2.8, mean(k15), .4, 'facecolor', colorVec(5,:));
hold on;
e5 = errorbar(2.8, mean(k15), std(k15), 'color', 'k', 'linewidth', 2);

h6 = bar(3.2, mean(k16), .4, 'facecolor', colorVec(6,:));
hold on;
e6 = errorbar(3.2, mean(k16), std(k16), 'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [1 2 3])
legend([h1 h2 h3 h4 h5 h6], {'mech + KSOM', 'mech + MM', ...
    'exp + KSOM', 'exp + MM', 'inc + KSOM', 'inc + MM'}, 'location', 'northwest');
set(gca, 'xticklabel', {'', '', '', '', '', ''});
xlim([.5 3.5]);
ylim([0 .22]);
grid on;
ylabel('stiffness ');
title('Oocyte stiffness under different conditions');

[p h] = ranksum(k11, k15)
[p h] = ranksum(k1Mech, k1Inc)
[p h] = ranksum(k1K, k1M)


%% 4. Plot mechanics over time

% group 1 = KSOM, group 2 = MM
% 1 = GV, 2 = M1, 3 = M2
colorVec = [0 0 .5; .6 .6 .8];
baseDir = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte';
groups = {'mechK', 'mechM', 'expK', 'expM', 'incK', 'incM'};
dates = {'10-8-15', '10-21-15'};
hourList = {'', '_8', '_16', '_24', '', '_8', '_16', '_24'};
groupAll = [];
dateAll = [];
embryoNumAll = [];
matStatus = [];

% 10-8-15 data
% which index of groups is embryo in?
groupAll = [groupAll 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2];
dateAll = [dateAll ones(1,20)];
embryoNumAll = [embryoNumAll 1:20];
matStatus = [2 3 3 3; NaN NaN NaN NaN; NaN NaN NaN NaN; 2 2 3 3; 2 2 3 3; ...
    2 2 3 3; 1 2 3 3; 2 2 3 3; 2 2 3 3; 1 2 2 2; ...
    1 1 2 2; NaN NaN NaN NaN; 2 2 3 3; 2 2 3 3; 2 2 2 2; ...
    1 2 2 2; 2 2 3 3; 2 2 2 2; 1 2 2 2; 1 2 3 3];

% 10-21-15
groupAll = [groupAll ones(1,10) 2*ones(1,20)];
dateAll = [dateAll 2*ones(1,30)];
embryoNumAll = [embryoNumAll 1:30];
matStatus = [matStatus; 1 1 1 1; 1 2 3 3; NaN NaN NaN NaN; 1 2 2 2; ...
    1 2 2 2; 1 1 1 1; 1 2 3 3; 1 2 2 2; 2 2 2 2; 1 2 3 3; 1 2 3 3; ...
    1 1 1 1; 1 3 3 3; 1 2 3 3; 1 2 3 3; 1 3 3 3; 1 1 1 1; 1 3 3 3; ...
    2 3 3 3; 2 3 3 3; 3 3 3 3; 1 1 1 1; NaN NaN NaN NaN; 2 2 3 3; ...
    2 3 3 3; 3 3 3 3; 2 3 3 3; 2 3 3 3; 2 2 3 3; 1 3 3 3];

% all mech params are stored in nx4 matrix
numEmbryos = length(groupAll);
k1list = zeros(numEmbryos,4);
n1list = zeros(numEmbryos,4);
k0list = zeros(numEmbryos,4);
taulist = zeros(numEmbryos,4);

% build vectors for each embryo
for i = 1:numEmbryos
    
    if ~isnan(groupAll(i))
        colorList(i,:) = colorVec(groupAll(i),:);
    end
    
    for j = 1:4
        extraString = hourList{4*(dateAll(i) - 1) + j};
        
        % bulk mech data
        dataDir = [baseDir '\' dates{dateAll(i)} ' analysis'];
        dataBulk = [dataDir, '\AutoMeasure\aspiration_data_', ...
            strrep(dates{dateAll(i)}, '-', '_'), ...
            '_E', num2str(embryoNumAll(i)), extraString, '.mat'];
        
        if ~exist(dataBulk, 'file')
            k1list(i,j) = NaN;
            n1list(i,j) = NaN;
            k0list(i,j) = NaN;
            taulist(i,j) = NaN;
        else
            load(dataBulk);
            k1list(i,j) = k1;
            n1list(i,j) = n1;
            k0list(i,j) = k0;
            taulist(i,j) = tau;
        end
    end
end

% k1list = k1list - repmat(k1list(:,1), 1, 4);

k1t1g1 = k1list(groupAll == 1 & ~isnan(k1list(:,1)'),1);
k1t2g1 = k1list(groupAll == 1 & ~isnan(k1list(:,2)'),2);
k1t3g1 = k1list(groupAll == 1 & ~isnan(k1list(:,3)'),3);
k1t4g1 = k1list(groupAll == 1 & ~isnan(k1list(:,4)'),4);
k1t1g2 = k1list(groupAll == 2 & ~isnan(k1list(:,1)'),1);
k1t2g2 = k1list(groupAll == 2 & ~isnan(k1list(:,2)'),2);
k1t3g2 = k1list(groupAll == 2 & ~isnan(k1list(:,3)'),3);
k1t4g2 = k1list(groupAll == 2 & ~isnan(k1list(:,4)'),4);

% bar graph over time
figure(1);
clf;
set(gca, 'fontsize', 14);
h1 = bar(.8, mean(k1t1g1), .4, 'facecolor', colorVec(1,:));
hold on;
e1 = errorbar(.8, mean(k1t1g1), std(k1t1g1), 'color', 'k', 'linewidth', 2);

h2 = bar(1.2, mean(k1t1g2), .4, 'facecolor', colorVec(2,:));
hold on;
e2 = errorbar(1.2, mean(k1t1g2), std(k1t1g2), 'color', 'k', 'linewidth', 2);

h3 = bar(1.8, mean(k1t2g1), .4, 'facecolor', colorVec(1,:));
hold on;
e3 = errorbar(1.8, mean(k1t2g1), std(k1t2g1), 'color', 'k', 'linewidth', 2);

h4 = bar(2.2, mean(k1t2g2), .4, 'facecolor', colorVec(2,:));
hold on;
e3 = errorbar(2.2, mean(k1t2g2), std(k1t2g2), 'color', 'k', 'linewidth', 2);

h5 = bar(2.8, mean(k1t3g1), .4, 'facecolor', colorVec(1,:));
hold on;
e5 = errorbar(2.8, mean(k1t3g1), std(k1t3g1), 'color', 'k', 'linewidth', 2);

h6 = bar(3.2, mean(k1t3g2), .4, 'facecolor', colorVec(2,:));
hold on;
e6 = errorbar(3.2, mean(k1t3g2), std(k1t3g2), 'color', 'k', 'linewidth', 2);

h7 = bar(3.8, mean(k1t4g1), .4, 'facecolor', colorVec(1,:));
hold on;
e7 = errorbar(3.8, mean(k1t4g1), std(k1t4g1), 'color', 'k', 'linewidth', 2);

h8 = bar(4.2, mean(k1t4g2), .4, 'facecolor', colorVec(2,:));
hold on;
e8 = errorbar(4.2, mean(k1t4g2), std(k1t4g2), 'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [1 2 3 4])
legend([h3 h4], {'KSOM', 'MM'}, 'location', 'northwest');
% legend([h1 h2 h3 h4 h5 h6 h7 h8], {'t = 0, KSOM', 't = 0, MM', ...
%     't = 8, KSOM', 't = 8, MM', 't = 16, KSOM', 't = 16, MM', ...
%     't = 24, KSOM', 't = 24, MM'}, 'location', 'northwest');
set(gca, 'xticklabel', {'t = 0', 't = 8', 't = 16', 't = 24'});
xlim([0.5 4.5]);
% ylim([-.03 .05]);
% ylim([.04 .16])
grid on;
ylabel('stiffness (k1)');
title('Oocyte stiffness during maturation');

[p h] = ranksum(k1t1g2,k1t4g2)
[p h] = ranksum(k1t1g1,k1t4g1)

hold on;
plot([.8 1.8 2.8 3.8] + .2, [mean(k1t1g1), mean(k1t2g1), mean(k1t3g1), mean(k1t4g1)], ...
    '--k', 'linewidth', 2);
plot([1.2 2.2 3.2 4.2] - .2, [mean(k1t1g2), mean(k1t2g2), mean(k1t3g2), mean(k1t4g2)], ...
    '--k', 'linewidth', 2);

% line plot over time
matColors = [.6 0 0; 0 0 .6; 0 .6 0]; % R = GV, B = M1, G = M2
figure(2);
subplot(1,2,1);
clf;
subplot(1,2,2);
clf;

% handles for legend
hList = zeros(2,3);

for i = 1:numEmbryos
    for j = 1:4
        
        % plot marker
        if ~isnan(matStatus(i,j))% && matStatus(i,4) == 3 
            
            subplot(1,2,groupAll(i));
            hold on;
%             hCurr = scatter(matStatus(i,j), k1list(i,j), 100, 'marker', 'o', ...
%                 'markerfacecolor', ...
%                 matColors(matStatus(i,j),:), 'markeredgecolor', 'none');
            hCurr = scatter(j, k1list(i,j), 100, 'marker', 'o', ...
                'markerfacecolor', ...
                matColors(matStatus(i,j),:), 'markeredgecolor', 'none');

            hList(groupAll(i), matStatus(i,j)) = hCurr;
            
        end
        
        % plot connecting line
        if (j < 4) && ~isnan(matStatus(i,j)) && ~isnan(matStatus(i,j+1)) % && matStatus(i,4) == 3
            subplot(1,2,groupAll(i));
            hold on;
%             line([matStatus(i,j) matStatus(i,j+1)], [k1list(i,j) k1list(i,j+1)], 'linestyle', ...
%                 '--', 'linewidth', 2, 'color', 'k');
            line([j j+1], [k1list(i,j) k1list(i,j+1)], 'linestyle', ...
                '--', 'linewidth', 2, 'color', 'k');
        end
        
    end
end

subplot(1,2,1);
set(gca, 'fontsize', 14);
title('oocyte maturation (KSOM)');
set(gca, 'xtick', [1 2 3 4]);
set(gca, 'xticklabel', {'0 hr', '8 hr', '16 hr', '24 hr'});
ylabel('stiffness (k1)');
% legend(hList(1,:), {'GV', 'MI', 'MII'}, 'location', 'northwest');
% ylim([-.04 .04]);
% ylim([.05 .12]);

subplot(1,2,2);
set(gca, 'fontsize', 14);
title('oocyte maturation (MM)');
set(gca, 'xtick', [1 2 3 4]);
set(gca, 'xticklabel', {'0 hr', '8 hr', '16 hr', '24 hr'});
ylabel('stiffness (k1)');
% legend(hList(2,:), {'GV', 'MI', 'MII'}, 'location', 'northwest');
% ylim([-.04 .04]);
% ylim([.05 .12]);

%% 5. Plot mechanics vs morphology at 24 hr time point

% morphology @ 24 hr time point
baseDir = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte';
dates = {'10-8-15', '10-21-15'};
groupAll = [];
mEnd = [];
dateAll = [];
embryoNumAll = [];

% 10-8-15 data
% which index of groups is embryo in?
% groupAll = [1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 5*ones(1,16) ...
%     6*ones(1,12) 4*ones(1,9) 3*ones(1,8)];
% mEnd = [3 1 NaN 3 3 3 3 3 3 2 2 NaN 3 3 2 2 3 2 2 3 2 2 3 3 2 ...
%     2 1 2 2 3 2 2 2 3 3 3 3 3 3 3 3 3 3 3 2 3 3 3 3 3 ...
%     3 2 2 2 1 2 2 3 3 3 2 2 2 2 2];
% dateAll = ones(1,65);
% embryoNumAll = 1:65;

% 10-21-15
groupAll = [groupAll ones(1,10) 2*ones(1,20) 3*ones(1,11) ...
    4*ones(1,10) 6*ones(1,10) 5*ones(1,10)];
mEnd = [mEnd 1 3 NaN 2 2 1 3 2 2 3 3 1 3 3 3 3 1 3 3 3 3 1 NaN 3 3 3 3 ...
    3 3 3 3 3 2 2 3 3 2 2 3 3 2 2 1 3 3 3 3 2 3 3 3 3 3 3 2 3 3 2 ...
    2 3 2 2 3 2 2 3 3 1 3 2 2];
dateAll = [dateAll 2*ones(1,71)];
embryoNumAll = [embryoNumAll 1:71];

colorVec = [.6 0 0; .6 .6 .8; 0 0 .6; .6 .8 .8; 0 .6 0; .6 .8 .6];

numEmbryos = length(mEnd);
k1list = zeros(1,numEmbryos);
n1list = zeros(1,numEmbryos);
k0list = zeros(1,numEmbryos);
taulist = zeros(1,numEmbryos);

% build vectors for each embryo
for i = 1:numEmbryos
    
    if ~isnan(mEnd(i))
        colorList(i,:) = colorVec(groupAll(i),:);
    end
    
    extraString = '';
    if groupAll(i) < 3
        extraString = '_24';
    end
    
    % bulk mech data
    dataDir = [baseDir '\' dates{dateAll(i)} ' analysis'];
    dataBulk = [dataDir, '\AutoMeasure\aspiration_data_', ...
        strrep(dates{dateAll(i)}, '-', '_'), ...
        '_E', num2str(embryoNumAll(i)), extraString, '.mat'];
    
    if ~exist(dataBulk, 'file')
        k1list(i) = NaN;
        n1list(i) = NaN;
        k0list(i) = NaN;
        taulist(i) = NaN;
    else
        load(dataBulk);
        k1list(i) = k1;
        n1list(i) = n1;
        k0list(i) = k0;
        taulist(i) = tau;
    end
end

% separate by morphology and culture conditions
kK1 = k1list(~isnan(k1list) & ~isnan(mEnd) & (mEnd == 1) & mod(groupAll,2) == 1);
kK2 = k1list(~isnan(k1list) & ~isnan(mEnd) & (mEnd == 2) & mod(groupAll,2) == 1);
kK3 = k1list(~isnan(k1list) & ~isnan(mEnd) & (mEnd == 3) & mod(groupAll,2) == 1);

kM1 = k1list(~isnan(k1list) & ~isnan(mEnd) & (mEnd == 1) & mod(groupAll,2) == 0);
kM2 = k1list(~isnan(k1list) & ~isnan(mEnd) & (mEnd == 2) & mod(groupAll,2) == 0);
kM3 = k1list(~isnan(k1list) & ~isnan(mEnd) & (mEnd == 3) & mod(groupAll,2) == 0);

xEst = .05:.0005:.15;

% make smoothed kernel density plot
[fk1 xk1] = ksdensity(kK1,xEst);
[fk2 xk2] = ksdensity(kK2,xEst);
[fk3 xk3] = ksdensity(kK3,xEst);
[fm1 xm1] = ksdensity(kM1,xEst);
[fm2 xm2] = ksdensity(kM2,xEst);
[fm3 xm3] = ksdensity(kM3,xEst);

figure(4);
clf;
subplot(2,1,1);
grid on;
set(gca, 'FontSize', 14);
title('Oocyte mechanics after maturation, KSOM');
xlabel('stiffness (k1)');
ylabel('histogram density');

% normalize hist densities to none stand out vertically too much
% fk1 = fk1 / max(fk1);
% fk2 = fk2 / max(fk2);
% fk3 = fk3 / max(fk3);
% fm1 = fm1 / max(fm1);
% fm2 = fm2 / max(fm2);
% fm3 = fm3 / max(fm3);

hold on;
hk1 = plot(xk1, fk1, 'Color', colorVec(1,:), 'linewidth', 2);
hk2 = plot(xk2, fk2, 'Color', colorVec(3,:), 'linewidth', 2);
hk3 = plot(xk3, fk3, 'Color', colorVec(5,:), 'linewidth', 2);

ak1 = area(xk1, fk1, 'EdgeColor', colorVec(1,:), 'FaceColor', colorVec(1,:));
ak1c = get(ak1, 'Children');
set(ak1c, 'FaceAlpha', .3);
ak2 = area(xk2, fk2, 'EdgeColor', colorVec(3,:), 'FaceColor', colorVec(3,:));
ak2c = get(ak2, 'Children');
set(ak2c, 'FaceAlpha', .3);
ak3 = area(xk3, fk3, 'EdgeColor', colorVec(5,:), 'FaceColor', colorVec(5,:));
ak3c = get(ak3, 'Children');
set(ak3c, 'FaceAlpha', .3);

legend([hk1 hk2 hk3], {'GV', 'MI', 'MII'}, 'box', 'off');
box on;
xlim([.05 .15]);

subplot(2,1,2);
grid on;
set(gca, 'FontSize', 14);
title('Oocyte mechanics after maturation, MM');
xlabel('stiffness (k1)');
ylabel('histogram density');

hold on;
hm1 = plot(xm1, fm1, 'Color', colorVec(1,:), 'linewidth', 2);
hm2 = plot(xm2, fm2, 'Color', colorVec(3,:), 'linewidth', 2);
hm3 = plot(xm3, fm3, 'Color', colorVec(5,:), 'linewidth', 2);

am1 = area(xm1, fm1, 'EdgeColor', colorVec(1,:), 'FaceColor', colorVec(1,:));
am1c = get(am1, 'Children');
set(am1c, 'FaceAlpha', .3);
am2 = area(xm2, fm2, 'EdgeColor', colorVec(3,:), 'FaceColor', colorVec(3,:));
am2c = get(am2, 'Children');
set(am2c, 'FaceAlpha', .3);
am3 = area(xm3, fm3, 'EdgeColor', colorVec(5,:), 'FaceColor', colorVec(5,:));
am3c = get(am3, 'Children');
set(am3c, 'FaceAlpha', .3);

legend([hm1 hm2 hm3], {'GV', 'MI', 'MII'}, 'box', 'off');
xlim([.05 .15]);


% make smoothed kernel density plot of KSOM vs MM only
xEst = .05:.0005:.15;
[fk xk] = ksdensity([kK1 kK2 kK3],xEst);
[fm xm] = ksdensity([kM1 kM2 kM3],xEst);

figure(5);
clf;
grid on;
set(gca, 'FontSize', 14);
title('Oocyte mechanics after maturation');
xlabel('stiffness (k1)');
ylabel('histogram density');

hold on;
hk = plot(xk, fk, 'Color', [0 0 .5], 'linewidth', 2);
hm = plot(xm, fm, 'Color', [.6 .6 .8], 'linewidth', 2);

ak = area(xk, fk, 'EdgeColor', [0 0 .5], 'FaceColor', [0 0 .5]);
akc = get(ak, 'Children');
set(akc, 'FaceAlpha', .3);
am = area(xm, fm, 'EdgeColor', [.6 .6 .8], 'FaceColor', [.6 .6 .8]);
amc = get(am, 'Children');
set(amc, 'FaceAlpha', .3);

legend([hk hm], {'KSOM', 'MM'}, 'box', 'off');
box on;


%% 6. For a single embryo, plot pressure vs equilibrium depth


load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\10-8-15 analysis\Zona\allZona.mat');
embryoToPlot = 10;

figure(6);
clf;
scatter((depthAll{embryoToPlot} - min(depthAll{embryoToPlot}))/20, ...
    (pressureAll{embryoToPlot} - min(pressureAll{embryoToPlot})), 100, ...
    'marker', '+', 'linewidth', 3, 'markerfacecolor', [0 0 .8]);
grid on;

hold on;
plot(ones(1,20), linspace(0,.5,20), '--k', 'linewidth', 3);
xlim([0 1.7]);
set(gca, 'fontsize', 14);


























