% track changes in mechanical properties over the course of
% oocyte maturation
% LZ 9-2-15


%% 1. Load all mechanical data

% oocyteNums should be [1 1 1 1 2 2 2 2 3 3 3 3] etc
% measHour should be [0 6 12 20 0 6 12 20 0 6 12 20] etc
% their length should be the number of measurements taken that day

numOocytesTotal = [];
numMeasTotal = [];
baseTime = [];
currDate = {};
oocyteNums = {};
measHour = {};

% trial 1: 8-27-15
measuredOocytes = repmat(1:8, 4, 1);
numOocytesTotal = [numOocytesTotal 8];
numMeasTotal = [numMeasTotal 32];
baseTime = [baseTime 8]; % # hrs after hCG injection 1st measurement done
currDate = {currDate{:}, '8-27-15'};
oocyteNums = {oocyteNums{:}, measuredOocytes(:)'}; % meas num for all oocytes
measHour = {measHour{:}, repmat([0 6 12 20],1,8)}; % meas time for all oocytes

filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\';
k0list = cell(1,length(currDate));
k1list = cell(1,length(currDate));
n0list = cell(1,length(currDate));
n1list = cell(1,length(currDate));
taulist = cell(1,length(currDate));
colorVec = zeros(length(currDate),3);

% C1 = 0; C2 = 6; C3 = 12; C4 = 20
colors = [.6 0 0; 0 .6 0; 0 0 .6; .6 0 .6];

for i = 1:length(currDate)
    
    currNumMeas = length(oocyteNums{i});
    
    % separate into cells by date
    k0list{i} = zeros(1,currNumMeas);
    k1list{i} = zeros(1,currNumMeas);
    n0list{i} = zeros(1,currNumMeas);
    n1list{i} = zeros(1,currNumMeas);
    taulist{i} = zeros(1,currNumMeas);
    
    for j = 1:currNumMeas
        
        oocyteString = num2str(oocyteNums{i}(j));
        
        if measHour{i}(j) == 0
            measString = '';
        else
            measString = ['_' num2str(measHour{i}(j))];
        end
        
        % try to find automated measurement first
        testPath = [filePath1, currDate{i}, ...
            ' analysis\AutoMeasure\aspiration_data_', ...
            strrep(currDate{i}, '-', '_'), '_E', oocyteString, ...
            measString, '.mat'];
        
        % then try manual measurement
        if ~exist(testPath, 'file')
            testPath = [filePath1, currDate{i}, ...
                ' analysis\aspiration_data_', ...
                strrep(currDate{i}, '-', '_'), '_E', oocyteString, ...
                measString, '.mat'];
        end
        
        % if it still doesn't exist, just put NaNs
        % otherwise, load parameters
        if ~exist(testPath, 'file')
            k0list{i}(j) = NaN;
            k1list{i}(j) = NaN;
            n0list{i}(j) = NaN;
            n1list{i}(j) = NaN;
            taulist{i}(j) = NaN;
        else
            load(testPath);
            k0list{i}(j) = k0;
            k1list{i}(j) = k1;
            n0list{i}(j) = n0;
            n1list{i}(j) = n1;
            taulist{i}(j) = tau;
        end
        
    end
end

%% 2. Gather parameter vectors to plot

% plot for each oocyte measured
plotHour = cell(1,sum(numOocytesTotal));
plotParam = cell(1,sum(numOocytesTotal));
currOocyte = 1;

for i = 1:length(currDate)
    currNumOocytes = numOocytesTotal(i);
    
    for j = 1:currNumOocytes
        
        currHourList = measHour{i};
        currOocyteNums = oocyteNums{i};
        currParamList = k1list{i};
        
        plotHour{currOocyte} = baseTime(i) + currHourList(currOocyteNums == j);
        plotParam{currOocyte} = currParamList(currOocyteNums == j);
        currOocyte = currOocyte + 1;
        
    end
end

%% 3. Plot

figure(1);
clf;
set(gca, 'fontsize', 14);
hold on;

for i = 1:sum(numOocytesTotal)
    
    plot(plotHour{i}, plotParam{i}, 'color', [0 0 1], ...
        'marker', 'o', 'linewidth', 2);
    text(7, plotParam{i}(1), num2str(i), 'fontsize', 14);
    text(29, plotParam{i}(end), num2str(i), 'fontsize', 14);

end

xlabel('time after hCG injection (hrs)');
ylabel('k1 parameter');
% ylim([.02 .12]);
xlim([6 30]);
title('Oocyte stiffness during in vitro maturation');









