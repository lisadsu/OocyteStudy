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
%                                    rest must have either 0 or 1
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


function [mOut, paramsOut, testIndList, enumList] = loadSVMdataClinical(inputMethod, ...
    nGroups, paramNumsToUse)


% ========================================
% MAKE SURE LIST OF EXPERIMENTS IS CORRECT
% ========================================
params = loadDataClinicalStudy();
numParticipants = params.numParticipants;
participantIDs = params.participantIDs;
numEmbryos = params.numEmbryos;
outcomeInfo = params.outcomeInfo;
zygoteM = params.zygoteM;
d3M = params.d3M;
blastM = params.blastM;
ICSI = params.ICSI;
PGD = params.PGD;
gender = params.gender;
totalNumEmbryos = sum(numEmbryos);
procDataPath = 'C:\Users\Livia\Desktop\IVF\Processed Data\Human\';

% init param vectors
mList = [];
k0list = [];
k1list = [];
n0list = [];
n1list = [];
taulist = [];
gList = [];
colorMat = [];
icsiList = [];
aList = [];
currE = 0;
zList = [];

ptsToPlot = ones(1,numParticipants);
aText = [];
cText = {};
dayThreeM1 = [];
dayThreeM2 = [];

% load in params by participant and embryo num
for i = 1:numParticipants
    if ptsToPlot(i)
        for j = 1:numEmbryos(i)
            
            currE = currE + 1;
            aText = [aText currE];
            cText = {cText{:}, ['P', num2str(i), '\_E', num2str(j)]};
            
            currDataPath = [procDataPath participantIDs{i} '\' ...
                participantIDs{i} '_E' num2str(j) '.mat'];
            
            % save embryo params and color
            if exist(currDataPath, 'file')% && (ICSI{i}(j) == 1)
                load(currDataPath);
                mList = [mList outcomeInfo{i}(j)];
                k0list = [k0list k0];
                k1list = [k1list k1];
                n0list = [n0list n0];
                n1list = [n1list n1];
                taulist = [taulist tau];
                gList = [gList gender{i}(j)];
                icsiList = [icsiList ICSI{i}(j)];
                zList = [zList zygoteM{i}(j)];
                
                aPad = padarray(A,[0,65-length(A)],'replicate', 'post');
                aList = [aList; aPad(1:65)];
                
                % get day 3 morphology data
                dayThreeM = d3M{i}{j};
                if numel(dayThreeM) == 0
                    dayThreeM1 = [dayThreeM1 NaN];
                    dayThreeM2 = [dayThreeM2 NaN];
                    k1list(end) = NaN;
                else
                    dayThreeM1 = [dayThreeM1 str2num(dayThreeM(1))];
                    dayThreeM2 = [dayThreeM2 getNumFromNumeral(dayThreeM(2:end))];
                end
                
                if isnan(outcomeInfo{i}(j))% || (ICSI{i}(j) == 1)
                    mList(end) = NaN;
                    k1list(end) = NaN;
                    colorMat = [colorMat; [NaN NaN NaN]];
                elseif outcomeInfo{i}(j)
                    colorMat = [colorMat; [0 .6 0]];
                else
                    colorMat = [colorMat; [0 0 .6]];
                end
                
            else
                colorMat = [colorMat; [NaN NaN NaN]];
                mList = [mList NaN];
                k0list = [k0list NaN];
                k1list = [k1list NaN];
                n0list = [n0list NaN];
                n1list = [n1list NaN];
                taulist = [taulist NaN];
                icsiList = [icsiList NaN];
                aList = [aList; NaN*zeros(1,65)];
                gList = [gList NaN];
                zList = [zList NaN];
                dayThreeM1 = [dayThreeM1 NaN];
                dayThreeM2 = [dayThreeM2 NaN];
            end
        end
    end
end


zList = zList + .3*randn(1,length(zList)); % add a little random noise to morphology
numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) & ...
        ~isnan(mList) & ~isnan(dayThreeM1);

if inputMethod == 2
    numToTest = length(mList(mList == 5));
    enumList = 1:numToTest;
    enumList = enumList(numsToPlot(end-numToTest+1:end));
else
    enumList = [];
end

% new names for variables
k = k1list(numsToPlot);
n = n1list(numsToPlot);
t = taulist(numsToPlot);
k2 = k0list(numsToPlot);
m1 = dayThreeM1(numsToPlot);
m2 = dayThreeM2(numsToPlot);
a = aList(numsToPlot,:);
z = zList(numsToPlot);
aMin = repmat(min(a, [], 1), size(a,1), 1);
aMax = repmat(max(a, [], 1), size(a,1), 1);

% normalize each parameter list from 0-1
kk = (k - min(k)) / (max(k) - min(k));
nn = (n - min(n)) / (max(n) - min(n));
tt = (t - min(t)) / (max(t) - min(t));
kk2 = (k2 - min(k2)) / (max(k2) - min(k2));
mm1 = (m1 - min(m1)) / (max(m1) - min(m1));
mm2 = (m2 - min(m2)) / (max(m2) - min(m2));
zz = (z - min(z)) / (max(z) - min(z));
aa = (a - aMin) ./ (aMax - aMin);

% Scale data first
% normalize all params 0-1
m = mList(numsToPlot);
paramsOut = [kk ; nn ; tt ; kk2; mm1; mm2; aa']'; %aa
paramsOut = paramsOut(:,paramNumsToUse);
paramsOut = [paramsOut];%, zz'];
mOut = m;

% ============================================================
% BEGINNING OF IF STATEMENT
% ============================================================

if inputMethod == 1
    
    % ignore embryos with unknown viability
    paramsOut = paramsOut(mOut ~= 2,:);
    mOut = mOut(mOut ~= 2);
    mOut(mOut > 0) = 1;
    
    % Create a n-fold cross-validation to compute classification error
    % m is ground truth data (+/-1)
    testIndList = crossvalind('Kfold', mOut, nGroups)';
    
elseif inputMethod == 2
    
    testIndList = ones(1,length(mOut));
    testIndList(mOut == 2) = 2;
    
end






