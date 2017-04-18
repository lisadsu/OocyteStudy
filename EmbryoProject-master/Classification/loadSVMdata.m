% load all data for running classification using SVM
% Inputs:
%   inputMethod: if inputMethod = 1, output training and test set by
%                                    randomly splitting data evenly
%
%                if inputMethod = 2, output training set as all data so
%                                    far, and test data as current
%                                    experiment (choose this option for
%                                    real-time classification)
%
%                                    the embryos to be classified must have
%                                    morphology values of 5, while all the
%                                    rest must have either 2 or 4
%
%   nGroups: Number of partitions to split data for validation
%            If inputMethod = 2, nGroups MUST equal 2
%            If inputMethod = 1, nGroups must be >= 2
%
% Outputs:
%   mOut = ground truth morphology data (contains test + train data)
%   paramsOut = param lists for test + train data
%   testIndList = indices indicating which measurement is in which group
%                 If inputMethod = 2, then testIndList = 1 is the training
%                       data and testIndList = 2 is the testing data
%                 If inputMethod = 1, then each value of testIndList
%                       represents a different partition of the data


function [mOut, paramsOut, testIndList, enumList] = loadSVMdata(inputMethod, ...
    nGroups, plotInput, filePath1, paramNumsToUse)

if plotInput
    figPlot = figure;
end

% ========================================
% MAKE SURE LIST OF EXPERIMENTS IS CORRECT
% ========================================

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
% whichToPlot(22) is 6-17-13 (fresh, used for cortical granule staining)
% whichToPlot(23) is 6-27-13 (fresh, used for cortical granule staining)

whichToPlot = [zeros(1,8) ones(1,5) 0 1 1 1 1 0 1 1 1 0];
plotAuto = [zeros(1,8) ones(1,15)];
pixelConv = [54*ones(1,10) 108*ones(1,length(whichToPlot)-10)];

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

dateList2 = cell(1,length(whichToPlot));
for j = 1:length(whichToPlot)
    dateU = dateList{j};
    dateUI = strfind(dateU, '-');
    dateU(dateUI) = '_';
    dateList2{j} = dateU;
end

% for making legend to scatter plot
colorMat = [];
k0list = [];
k1list = [];
n0list = [];
taulist = [];
elonglist = [];
n1list = [];
morphologyOut = [];
timeLapseOut = [];

for jNum = 1:length(whichToPlot)
    if whichToPlot(jNum)
        
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
        
        colorMat = [colorMat ; cCurr];
        
        % add on to morphologyOut
        % group 1 and 2 (no blast and poor quality blast) into nonviable class
        % group 3 and 4 (medium and good quality blast) into viable class
        morphologyAdd = morphology{jNum};
        morphologyAdd(morphologyAdd == 3) = 4;
        morphologyAdd(morphologyAdd < 3) = 2;

        % add to morphologyOut
        morphologyOut = [morphologyOut morphologyAdd];        
        
        % add to timeLapseOut
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
                
                load([filePath1 filePath2 '_E', embryoString, '.mat']);
                
                aspiration_depth = aspiration_depth * 40 * 10^-6 / pixelConv(jNum); % convert from pixels to meters
                k0list = [k0list k0];
                k1list = [k1list k1];
                n0list = [n0list n0];
                taulist = [taulist tau];
                n1list = [n1list n1];
                elonglist = [elonglist F0/(k0 + k1)];
                
                currColor = cCurr(iNum,:);
                xdata = t;
                ydata = aspiration_depth;
                
                if plotInput
                    figure(figPlot);
                    plot(xdata, ydata, 'ob', 'Color', currColor);
                    hold on;
                    plot(xfine - min(xfine), yfit, 'Color', currColor);
                end
                
            else
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

if plotInput
    xlim([0 .45]);
    set(gca, 'FontSize', 14);
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
    title('Aspiration Depth of Mouse Embryos');
end

% put all variables and labels together in a list
% get rid of NaNs (morphologyOut has some NaNs, and so does k1list, etc
% but in different locations

% if not using any time lapse params, can use all mechanical data
if isempty(intersect(paramNumsToUse,[5 6 7]))
    numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) & ...
        ~isnan(morphologyOut);
else
    numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) & ...
        ~isnan(morphologyOut) & ~isnan(timeLapseOut(:,1))';
end

if inputMethod == 2
    numToTest = length(morphologyOut(morphologyOut == 5));
    enumList = 1:numToTest;
    enumList = enumList(numsToPlot(morphologyOut == 5));
else
    enumList = [];
end

% new names for variables
k = k1list(numsToPlot);
n = n1list(numsToPlot);
t = taulist(numsToPlot);
k2 = k0list(numsToPlot);

n = log(n);
k2 = log(k2);

% normalize each parameter list from 0-1
kk = (k - min(k)) / (max(k) - min(k));
nn = (n - min(n)) / (max(n) - min(n));
tt = (t - min(t)) / (max(t) - min(t));
kk2 = (k2 - min(k2)) / (max(k2) - min(k2));

% Scale data first
% normalize all params 0-1
m = morphologyOut(numsToPlot);
m(m == 2) = 0;
m(m == 4) = 1;

% [k1, n1, tau, k0, TL1, TL2, TL3]
paramsOut = [kk ; nn ; tt ; kk2]';

% process time lapse parameters and add to paramsOut
if ~isempty(intersect(paramNumsToUse,[5 6 7]))
    
    timeLapseOut = timeLapseOut(numsToPlot,:);
    
    timeLapseOut(:,1) = (timeLapseOut(:,1) - min(timeLapseOut(:,1))) / ...
        (max(timeLapseOut(:,1)) - min(timeLapseOut(:,1)));
    timeLapseOut(:,2) = (timeLapseOut(:,2) - min(timeLapseOut(:,2))) / ...
        (max(timeLapseOut(:,2)) - min(timeLapseOut(:,2)));
    timeLapseOut(:,3) = (timeLapseOut(:,3) - min(timeLapseOut(:,3))) / ...
        (max(timeLapseOut(:,3)) - min(timeLapseOut(:,3)));
    
    paramsOut = [paramsOut timeLapseOut];
end

paramsOut = paramsOut(:,paramNumsToUse);
mOut = m;

% ============================================================size(para
% BEGINNING OF IF STATEMENT
% ============================================================

if inputMethod == 1
    
    % ignore embryos with unknown viability
    paramsOut = paramsOut(m < 2,:);
    m = m(m < 2);
    
    % Create a n-fold cross-validation to compute classification error
    % m is ground truth data (+/-1)
    testIndList = crossvalind('Kfold', m, nGroups)';
    
elseif inputMethod == 2
    
    testIndList = ones(1,length(mOut));
    testIndList(mOut > 1) = 2;
    
end

end


