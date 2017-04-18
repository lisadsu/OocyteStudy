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


function [mOut, paramsOut, testIndList, enumList] = loadSVMdataOocyte(inputMethod, ...
    nGroups, paramNumsToUse)

groundTruth = 'fertList';

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
currE = 0;
colorList = [[.2 .2 .6];[.8 .2 .8];[.6 .6 .6]];

exptsToPlot = [ones(1,12)];
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
            if exist(currDataPath, 'file') && (morphologyInfo{i}(j) > -1) &&  measHour{i}(j) == 14
                
                load(currDataPath);
                
                fertList = [fertList fertInfo{i}(j)];
                blastList = [blastList blastForm{i}(j)];
                hatchList = [hatchList hatchInfo{i}(j)];
                k0list = [k0list k0];
                k1list = [k1list (k1+k1ScaleFactor(i))];
                n0list = [n0list n0];
                n1list = [n1list n1];
                taulist = [taulist tau];
                oocyteM = [oocyteM morphologyInfo{i}(j)];
                matEnv = [matEnv maturationEnv{i}(j)];
                aText = [aText currE];
                measHourList = [measHourList measHour{i}(j)];
                cText = {cText{:}, [currDate, '\_E', num2str(oocyteNums{i}(j)) params.fileNameApp{i}{j}]};
                
                aPad = padarray(A,[0,65-length(A)], 'replicate', 'post');
                aList = [aList; aPad(1:65)];
                
                
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

% put all variables and labels together in a list
% get rid of NaNs (morphologyOut has some NaNs, and so does k1list, etc
% but in different locations
groundTruth = eval(groundTruth);
numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(oocyteM) & ~isnan(groundTruth);

% new names for variables
k = k1list(numsToPlot);
n = n1list(numsToPlot);
t = taulist(numsToPlot);
k2 = k0list(numsToPlot);
m1 = oocyteM(numsToPlot);


% normalize each parameter list from 0-1
kk = (k - min(k)) / (max(k) - min(k));
nn = (n - min(n)) / (max(n) - min(n));
tt = (t - min(t)) / (max(t) - min(t));
kk2 = (k2 - min(k2)) / (max(k2) - min(k2));
mm1 = m1/2; %(m1 - min(m1)) / (max(m1) - min(m1));

% Scale data first
% normalize all params 0-1
m = groundTruth(numsToPlot);
paramsOut = [kk ; nn ; tt ; kk2; mm1]'; %aa
paramsOut = paramsOut(:,paramNumsToUse);
mOut = m;

% ============================================================
% BEGINNING OF IF STATEMENT
% ============================================================


if inputMethod == 1
    
    enumList = [];
    
    % ignore embryos with unknown viability
    paramsOut = paramsOut(mOut ~= 2,:);
    mOut = mOut(mOut ~= 2);
    mOut(mOut > 0) = 1;
    
    % Create a n-fold cross-validation to compute classification error
    % m is ground truth data (+/-1)
    testIndList = crossvalind('Kfold', mOut, nGroups)';
    
elseif inputMethod == 2
    
    numToTest = length(groundTruth(groundTruth == 2));
    enumList = 1:numToTest;
    enumList = enumList(numsToPlot(end-numToTest+1:end));
    
    testIndList = ones(1,length(mOut));
    testIndList(mOut == 2) = 2;
    
end



