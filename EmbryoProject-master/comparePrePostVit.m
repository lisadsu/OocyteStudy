% compare embryo mechanics pre and post-vitrification
% Livia Zarnescu
% 7-13-15

%% 1. Build parameter vectors for pre and post vit

clear all;
numEmbryos = 0;
currDate = {};
embryoNums = {};
exptCount = [];
measHour = {};
measMins = {};
measGroup = {};
% embryoNames = {};

% trial 1: 7-9-15
% embryos to exclude b/c they died: none
numEmbryos = numEmbryos + 5;
currDate = {currDate{:}, '7-9-15'};
embryoNums = {embryoNums{:}, [1 3 4 5 6]};
measHour = {measHour{:}, [0 1 2], [2 3 4], [2 3 4], [3 4 6], [3 4 6]};
measMins = {measMins{:}, {'44', '23', '08'}, {'48', '34', '33'}, {'54', '35', '34'}, ...
    {'39', '35', '10'}, {'40', '37', '11'}};
measGroup = {measGroup{:}, [4 4 4 4 4]};

% trial 2: 7-14-15
% embryos to exclude b/c they died: 4
numEmbryos = numEmbryos + 6;
currDate = {currDate{:}, '7-14-15'};
embryoNums = {embryoNums{:}, [1 2 3 5 8 10]};
measHour = {measHour{:}, [3 4 6], [3 4 6], [3 5 6], ...
    [4 5 6], [4 5 6], [4 5 6]};
measMins = {measMins{:}, {'41', '57', '11'}, {'42', '59', '13'}, {'43', '00', '14'}, ...
    {'08', '03', '16'}, {'09', '05', '18'}, {'10', '07', '20'}};
measGroup = {measGroup{:}, [4 4 4 4 4 4]};


% trial 3: 7-17-15
% embryos to exclude b/c they died: 12,13,25,28
numEmbryos = numEmbryos + 19;
currDate = {currDate{:}, '7-17-15'};
embryoNums = {embryoNums{:}, [1 2 3 4 5 6 9 10 14 17 18 19 20 21 22 26 27 29 30]};
measHour = {measHour{:}, [5], [5], [5], [6], [6], [6], [5], [5], ...
    [6], [2 4 5], [2 4 5], [2 4 5], [4 5 6], [4 5 6], [4 5 6], ...
    [2 4 5], [3 4 5], [4 5 6], [4 5 6]};
measMins = {measMins{:}, {'20'}, {'22'}, {'23'}, {'30'}, {'32'}, {'33'}, ...
    {'25'}, {'26'}, {'37'}, {'54', '30', '28'}, ...
    {'56', '32', '28'}, {'57', '33', '30'}, {'16', '47', '37'}, ...
    {'17', '46', '39'}, {'18', '45', '40'}, ...
    {'59', '35', '33'}, {'01', '37', '31'}, ...
    {'21', '49', '41'}, {'22', '50', '43'}};
measGroup = {measGroup{:}, [1 1 1 1 1 1 3 3 3 2 2 2 2 2 2 2 4 4 4 4]};

% trial 4: 7-20-15
% embryos to exclude b/c they died: 
numEmbryos = numEmbryos + 17;
currDate = {currDate{:}, '7-20-15'};
embryoNums = {embryoNums{:}, [1 2 3 4 5 6 7 9 13 14 15 16 17 18 20 21 23]};
measHour = {measHour{:}, [6], [6], [6], [6], [6], [6], [6], [6], ...
    [3 5 6], [3 5 6], [3 5 6], [6 7], [5 6 7], [5 6 7], ...
    [5 6 7], [3 5 6], [3 5 6]};
measMins = {measMins{:}, {'16'}, {'18'}, {'19'}, {'20'}, {'21'}, {'21'}, ...
    {'23'}, {'24'}, {'54', '17', '34'}, {'55', '15', '35'}, ...
    {'56', '14', '37'}, {'38', '14'}, {'01', '39', '16'}, ...
    {'02', '40', '17'}, {'04', '41', '18'}, {'57', '13', '42'}, ...
    {'58', '12', '44'}};
measGroup = {measGroup{:}, [1 1 1 1 1 1 3 3 2 2 2 2 2 2 4 4 4]};

% trial 5: 7-30-15
% embryos to exclude b/c they died: 
numEmbryos = numEmbryos + 16;
currDate = {currDate{:}, '7-30-15'};
embryoNums = {embryoNums{:}, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]};
measHour = {measHour{:}, [2], [2], [2], [2], [2], [3], [3], [3], ...
    [6], [6], [6], [6], [6], [6], [6], [6]};
measMins = {measMins{:}, {'53'}, {'55'}, {'56'}, {'58'}, {'59'}, {'00'}, ...
    {'01'}, {'02'}, {'24'}, {'26'}, {'27'}, {'28'}, {'29'}, {'30'}, ...
    {'32'}, {'33'}};
measGroup = {measGroup{:}, [1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3]};

% trial 6: 7-31-15
% embryos to exclude b/c they died: 11,12,14
numEmbryos = numEmbryos + 13;
currDate = {currDate{:}, '7-31-15'};
embryoNums = {embryoNums{:}, [1 2 3 4 5 6 7 8 9 10 13 15 16]};
measHour = {measHour{:}, [0 1 3], [0 1 3], [0 1 3], [0 1 3], ...
    [0 1 3], [0 1 3], [0 1 3], [0 1 3], ...
    [1 2 4], [1 2 4], [1 2 4], [1 2 4], [1 2 4]};
measMins = {measMins{:}, {'21', '36', '06'}, {'22', '37', '04'}, ...
    {'23', '38', '03'}, {'24', '39', '01'}, {'25', '41', '10'}, ...
    {'27', '42', '09'}, {'28', '43', '08'}, {'29', '44', '07'}, ...
    {'23', '49', '22'}, {'24', '50', '20'}, {'25', '52', '19'}, ...
    {'27', '53', '18'}, {'28', '54', '17'}};
measGroup = {measGroup{:}, [2 2 2 2 2 2 2 2 4 4 4 4 4]};


filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse embryo\';
vitStatus = {'preVit', 'postVit'};

k0list = cell(1,numEmbryos);
k1list = cell(1,numEmbryos);
n0list = cell(1,numEmbryos);
n1list = cell(1,numEmbryos);
taulist = cell(1,numEmbryos);
groupList = zeros(1,numEmbryos);
colorVec = zeros(numEmbryos,3);
embryoNumsAll = zeros(1, numEmbryos);
currEmbryo = 0;

% G1 = red; G2 = green; G3 = blue; G4 = purple
colors = [.6 0 0; 0 .6 0; 0 0 .6; .6 0 .6];


for n = 1:length(currDate)
    for j = 1:length(embryoNums{n})
      
        currEmbryo = currEmbryo + 1;
        colorVec(currEmbryo,:) = colors(measGroup{n}(j), :);
        groupList(currEmbryo) = measGroup{n}(j);
        embryoNumsAll(currEmbryo) = embryoNums{n}(j);
        
        for i = 1:2
            
            % load each file for each embryo
            embryoString = num2str(embryoNums{n}(j));
            numMeas = length(measHour{currEmbryo});
            
            if i == 1
                
                k0list{currEmbryo} = zeros(1,numMeas+1);
                k1list{currEmbryo} = zeros(1,numMeas+1);
                n0list{currEmbryo} = zeros(1,numMeas+1);
                n1list{currEmbryo} = zeros(1,numMeas+1);
                taulist{currEmbryo} = zeros(1,numMeas+1);
                
                % pre-vit
                testPath = [filePath1, currDate{n}, ...
                    ' analysis\preVit\AutoMeasure\aspiration_data_', ...
                    strrep(currDate{n}, '-', '_'), '_E', embryoString, '.mat'];
                
                if ~exist(testPath, 'file')
                    testPath = [filePath1, currDate{n}, ...
                        ' analysis\preVit\aspiration_data_', ...
                        strrep(currDate{n}, '-', '_'), '_E', embryoString, '.mat'];
                end
                
                load(testPath);
                k0list{currEmbryo}(1) = k0;
                k1list{currEmbryo}(1) = k1;
                n0list{currEmbryo}(1) = n0;
                n1list{currEmbryo}(1) = n1;
                taulist{currEmbryo}(1) = tau;
                
            elseif i == 2
                
                for k = 1:numMeas
                    
                    embryoTime = [num2str(measHour{currEmbryo}(k)) measMins{currEmbryo}{k}];
                    
                    % post-vit
                    testPath = [filePath1, currDate{n}, ...
                        ' analysis\postVit\AutoMeasure\aspiration_data_', ...
                        strrep(currDate{n}, '-', '_'), '_E', embryoString, '_', ...
                        embryoTime, '.mat'];
                    
                    if ~exist(testPath, 'file')
                        testPath = [filePath1, currDate{n}, ...
                            ' analysis\postVit\aspiration_data_', ...
                            strrep(currDate{n}, '-', '_'), '_E', embryoString, '_', ...
                            embryoTime, '.mat'];
                    end
                    
                    load(testPath);
                    k0list{currEmbryo}(k+1) = k0;
                    k1list{currEmbryo}(k+1) = k1;
                    n0list{currEmbryo}(k+1) = n0;
                    n1list{currEmbryo}(k+1) = n1;
                    taulist{currEmbryo}(k+1) = tau;
                    
                end
                
            end
        end
    end
end



%% 2. Build time list and plot vectors

% first, make list of measurement times
numEmbryosTotal = length(k0list);
timeList = cell(1,numEmbryosTotal);

% assign plot times on x axis
for i = 1:numEmbryosTotal
    
    if groupList(i) == 1 || groupList(i) == 3
        % one meas only at 0 and 3 hrs
        timeList{i} = [-5, 180];
        
    elseif groupList(i) == 2 || groupList(i) == 4
        
        % meas at 0 (pre-thaw), and (.5, 1.5, 3) post-thaw
        timeList{i} = zeros(1,length(measHour{i}) + 1);
        timeList{i}(1) = -5;
        
        for j = 1:length(measHour{i})
            timeList{i}(j+1) = measHour{i}(j)*60 + str2double(measMins{i}{j});
        end
        
        startTime = timeList{i}(2);
        
        for j = 1:length(measHour{i})
            timeList{i}(j+1) = timeList{i}(j+1) - startTime + 20;
        end
        
        
    end
    
end

%% 3. Plot parameters over time

figure(1);
clf;
firstOfGroup = [0 0 0 0];
legendNames = {'fresh, meas [0,3] hr', 'fresh, meas [0,.5,1.5,3] hr', ...
    'freeze/thaw, meas [0,3] hr', 'freeze/thaw, meas [0,.5,1.5,3] hr'};

for i = 1:numEmbryosTotal
    h(i) = plot(timeList{i}, k1list{i}, '-o', 'linewidth', 2, 'color', colorVec(i,:));
    hold on;
    
    if firstOfGroup(groupList(i)) == 0
        firstOfGroup(groupList(i)) = h(i);
    end
end


plot([0 0], [.05 .25], '--', 'color', [0 0 0], 'linewidth', 2);

set(gca, 'fontsize', 14);
legend(firstOfGroup(firstOfGroup ~= 0), legendNames(firstOfGroup ~= 0), ...
    'location', 'southeast');
xlabel('time post thaw (minutes)');
ylabel('k1 parameter (stiffness)');
title('Embryo mechanics after thawing');
xlim([-5 180]);
ylim([.05 .23]);
grid on;


%% 4. Plot changes due to freeze/thaw and repeated measurements

% Figure 2: average increase in stiffness over 3 hrs
preThaw = [];
postThaw = [];

for i = 1:numEmbryosTotal
    preThaw = [preThaw k1list{i}(1)];
    postThaw = [postThaw k1list{i}(end)];
end

thawDiff = postThaw - preThaw;
meanDiffGroup1 = mean(thawDiff(groupList == 1));
% thawDiff = thawDiff - meanDiffGroup1;

figure(2);
clf;

% Group 1
h11 = bar(.9, mean(thawDiff(groupList == 1)), .7, 'facecolor', [.8 .8 .8]);
hold on;
e11 = terrorbar(.9, mean(thawDiff(groupList == 1)),...
    std(thawDiff(groupList == 1)), .1);
set(e11, 'color', 'k', 'linewidth', 2);

% h12 = bar(1.1, mean(postThaw(groupList == 1)), .2, 'facecolor', [.3 .3 .3]);
% hold on;
% e12 = errorbar(1.1, mean(postThaw(groupList == 1)),...
%     std(postThaw(groupList == 1)),'color', 'k', 'linewidth', 2);
% setErrorBar(e12, 1.1, .03);

% Group 2
h21 = bar(1.9, mean(thawDiff(groupList == 2)), .7, 'facecolor', [.8 .8 .8]);
hold on;
e21 = terrorbar(1.9, mean(thawDiff(groupList == 2)),...
    std(thawDiff(groupList == 2)), .1);
set(e21, 'color', 'k', 'linewidth', 2);

% h22 = bar(2.1, mean(postThaw(groupList == 2)), .2, 'facecolor', [.3 .3 .3]);
% hold on;
% e22 = errorbar(2.1, mean(postThaw(groupList == 2)),...
%     std(postThaw(groupList == 2)),'color', 'k', 'linewidth', 2);
% setErrorBar(e22, 2.1, .03);

% Group 3
h31 = bar(2.9, mean(thawDiff(groupList == 3)), .7, 'facecolor', [.8 .8 .8]);
hold on;
e31 = terrorbar(2.9, mean(thawDiff(groupList == 3)),...
    std(thawDiff(groupList == 3)), .1);
set(e31, 'color', 'k', 'linewidth', 2);

% h32 = bar(3.1, mean(postThaw(groupList == 3)), .2, 'facecolor', [.3 .3 .3]);
% hold on;
% e32 = errorbar(3.1, mean(postThaw(groupList == 3)),...
%     std(postThaw(groupList == 3)),'color', 'k', 'linewidth', 2);
% setErrorBar(e32, 3.1, .03);

% Group 4
h41 = bar(3.9, mean(thawDiff(groupList == 4)), .7, 'facecolor', [.8 .8 .8]);
hold on;
e41 = terrorbar(3.9, mean(thawDiff(groupList == 4)),...
    std(thawDiff(groupList == 4)), .1);
set(e41, 'color', 'k', 'linewidth', 2);

% h42 = bar(4.1, mean(postThaw(groupList == 4)), .2, 'facecolor', [.3 .3 .3]);
% hold on;
% e42 = errorbar(4.1, mean(postThaw(groupList == 4)),...
%     std(postThaw(groupList == 4)),'color', 'k', 'linewidth', 2);
% setErrorBar(e42, 4.1, .03);


set(gca, 'xtick', [1 2 3 4])
set(gca, 'fontsize', 14);
% ylim([.05 .2]);
set(gca, 'xticklabel', {'Group1', 'Group2', 'Group3', 'Group4'});
ylabel('k1 parameter (stiffness)');
title('Change in stiffness after 3 hours');
xlim([0.5 4.5]);

[h p] = ttest2(thawDiff(groupList == 1), thawDiff(groupList == 3))
[h p] = ttest2(thawDiff(groupList == 2), thawDiff(groupList == 4))

%% Plot changes over time due to vit


% Figure 3: stiffness over time
timeData = [];
timeGroup = [];


for i = 1:numEmbryosTotal
    if (length(k1list{i}) == 4)
        timeData = [timeData; [0 k1list{i}(2:4) - k1list{i}(1)]];
        timeGroup = [timeGroup groupList(i)];
    end
end

figure(3);
clf;

% % time 1
% h12 = bar(0.9, mean(timeData(timeGroup == 2,1)), .2, 'facecolor', [.3 .3 .3]);
% hold on;
% e12 = errorbar(0.9, mean(timeData(timeGroup == 2,1)),...
%     std(timeData(timeGroup == 2,1)),'color', 'k', 'linewidth', 2);
% setErrorBar(e12, 0.9, .03);
% 
% h14 = bar(1.1, mean(timeData(timeGroup == 4,1)), .2, 'facecolor', [.6 .6 .6]);
% hold on;
% e14 = errorbar(1.1, mean(timeData(timeGroup == 4,1)),...
%     std(timeData(timeGroup == 4,1)),'color', 'k', 'linewidth', 2);
% setErrorBar(e14, 1.1, .03);
% 
% h11 = bar(0.7, mean(preThaw(groupList == 1)), .2, 'facecolor', [.1 .1 .1]);
% hold on;
% e11 = errorbar(0.7, mean(preThaw(groupList == 1)),...
%     std(preThaw(groupList == 1)),'color', 'k', 'linewidth', 2);
% setErrorBar(e11, 0.7, .03);
% 
% h13 = bar(1.3, mean(preThaw(groupList == 3)), .2, 'facecolor', [.9 .9 .9]);
% hold on;
% e13 = errorbar(1.3, mean(preThaw(groupList == 3)),...
%     std(preThaw(groupList == 3)),'color', 'k', 'linewidth', 2);
% setErrorBar(e13, 1.3, .03);

% time 2
h22 = bar(1.9, mean(timeData(timeGroup == 2,2)), .2, 'facecolor', [.3 .3 .3]);
hold on;
e22 = terrorbar(1.9, mean(timeData(timeGroup == 2,2)),...
    std(timeData(timeGroup == 2,2)), .1);
set(e22, 'color', 'k', 'linewidth', 2);

h24 = bar(2.1, mean(timeData(timeGroup == 4,2)), .2, 'facecolor', [.6 .6 .6]);
hold on;
e24 = terrorbar(2.1, mean(timeData(timeGroup == 4,2)),...
    std(timeData(timeGroup == 4,2)), .1);
set(e24, 'color', 'k', 'linewidth', 2);

% time 3
h32 = bar(2.9, mean(timeData(timeGroup == 2,3)), .2, 'facecolor', [.3 .3 .3]);
hold on;
e32 = terrorbar(2.9, mean(timeData(timeGroup == 2,3)),...
    std(timeData(timeGroup == 2,3)), .1);
set(e32, 'color', 'k', 'linewidth', 2);

h34 = bar(3.1, mean(timeData(timeGroup == 4,3)), .2, 'facecolor', [.6 .6 .6]);
hold on;
e34 = terrorbar(3.1, mean(timeData(timeGroup == 4,3)),...
    std(timeData(timeGroup == 4,3)), .1);
set(e34, 'color', 'k', 'linewidth', 2);

% time 4
h42 = bar(3.9, mean(timeData(timeGroup == 2,4)), .2, 'facecolor', [.3 .3 .3]);
hold on;
e42 = terrorbar(3.9, mean(timeData(timeGroup == 2,4)),...
    std(timeData(timeGroup == 2,4)), .1);
set(e42, 'color', 'k', 'linewidth', 2);

h44 = bar(4.1, mean(timeData(timeGroup == 4,4)), .2, 'facecolor', [.6 .6 .6]);
hold on;
e44 = terrorbar(4.1, mean(timeData(timeGroup == 4,4)),...
    std(timeData(timeGroup == 4,4)), .1);
set(e44, 'color', 'k', 'linewidth', 2);

h41 = bar(3.7, mean(thawDiff(groupList == 1)), .2, 'facecolor', [.1 .1 .1]);
hold on;
e41 = terrorbar(3.7, mean(thawDiff(groupList == 1)),...
    std(thawDiff(groupList == 1)), .1);
set(e41, 'color', 'k', 'linewidth', 2);

h43 = bar(4.3, mean(thawDiff(groupList == 3)), .2, 'facecolor', [.9 .9 .9]);
hold on;
e43 = terrorbar(4.3, mean(thawDiff(groupList == 3)),...
    std(thawDiff(groupList == 3)), .1);
set(e43, 'color', 'k', 'linewidth', 2);

plot([1 5], mean(thawDiff(groupList == 1))*ones(1,2), '--k', 'linewidth', 2);

set(gca, 'xtick', [2 3 4])
set(gca, 'fontsize', 14);
ylim([-.02 .10]);
legend([h41 h42 h44 h43], {'Fresh, meas 1x', ...
    'Fresh, meas 3x', 'Freeze/thaw, meas 3x', 'Freeze/thaw, meas 1x'}, ...
    'location', 'northwest');
set(gca, 'xticklabel', {'30 min', '120 min', '180 min'});
ylabel('change in k1 parameter (stiffness)');
xlabel('time post-thaw');
title('Change in stiffness after vitrification');
xlim([1.5 4.5]);
grid on;

[h p] = ttest2(timeData(timeGroup == 4,2),timeData(timeGroup == 2,2))


%% Make Supplementary Figure for paper

figure(4);
clf;

% 30 min, control group
h1c = bar(1.9, mean(timeData(timeGroup == 2,2)), .2, 'facecolor', [.1 .1 .1]);
hold on;
e1c = terrorbar(1.9, mean(timeData(timeGroup == 2,2)),...
    std(timeData(timeGroup == 2,2)), .1);
set(e1c, 'color', 'k', 'linewidth', 2);


% 30 min, meas group
h1m = bar(2.1, mean(timeData(timeGroup == 4,2)), .2, 'facecolor', [.7 .7 .7]);
hold on;
e1m = terrorbar(2.1, mean(timeData(timeGroup == 4,2)),...
    std(timeData(timeGroup == 4,2)), .1);
set(e1m, 'color', 'k', 'linewidth', 2);

% 3 hr, control group
h2c = bar(2.9, mean(thawDiff(groupList == 1)), .2, 'facecolor', [.1 .1 .1]);
hold on;
e2c = terrorbar(2.9, mean(thawDiff(groupList == 1)),...
    std(thawDiff(groupList == 1)), .1);
set(e2c, 'color', 'k', 'linewidth', 2);

% 3 hr, meas group
h2m = bar(3.1, mean(thawDiff(groupList == 3)), .2, 'facecolor', [.7 .7 .7]);
hold on;
e2m = terrorbar(3.1, mean(thawDiff(groupList == 3)),...
    std(thawDiff(groupList == 3)), .1);
set(e2m, 'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [2 3])
set(gca, 'fontsize', 14);
ylim([-.01 .08]);
legend([h1c h1m], {'Fresh (control)', 'Frozen/thawed'}, ...
    'location', 'northwest');
set(gca, 'xticklabel', {'30 min', '3 hr'});
ylabel('change in k_1 parameter from t_0');
xlabel('time post-thaw');
title('Change in stiffness after vitrification');
xlim([1.5 3.5]);
grid on;

[h p] = ttest2(timeData(timeGroup == 4,2),timeData(timeGroup == 2,2))
[h p] = ttest2(thawDiff(groupList == 1), thawDiff(groupList == 3))









