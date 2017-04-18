%
% load in either mouse embryo or oocyte params (in order: k1, n1, tau, k0)
%
% input params: type = 'mouse oocyte' or 'mouse embryo'
%               whichToPlot = vector of 0s and 1s selecting which
%                             experiments to include in output params
%
% output params: paramsOut = nx4 matrix, with columns: k1, n1, tau, k0
%                morphologyOut = nx1 vector of viability values


function [paramsOut, morphologyOut] = getMouseParams(type, whichToPlot)


% load params based on type
if strcmp(type, 'mouse embryo')
    
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

    plotAuto = [zeros(1,8) ones(1, length(whichToPlot)-8)];
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
    
    for jNum = 1:length(whichToPlot)
        if whichToPlot(jNum)
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
            cCurr(morphology{jNum} == 6, :) = repmat([0 .6 .6], s(1), 1);
            
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
                        
        end
    end

    numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) & ...
        ~isnan(mList);% & ~isnan(timeLapseOut(:,1))';
    mN = mList(numsToPlot);
    k1N = k1list(numsToPlot);
    n1N = n1list(numsToPlot);
    tN = taulist(numsToPlot);
    k0N = k0list(numsToPlot);
    
    paramsOut = [k1N', n1N', tN', k0N'];
    morphologyOut = mN';
    
elseif strcmp(type, 'mouse oocyte')
    
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
    
    plotAuto = [ones(1,length(whichToPlot))];
    pixelConv = [108*ones(1,length(whichToPlot))];

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
                        
        end
    end
    
    
    numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) & ~isnan(mList);
    mN = mList(numsToPlot);
    k1N = k1list(numsToPlot);
    n1N = n1list(numsToPlot);
    tN = taulist(numsToPlot);
    k0N = k0list(numsToPlot);
    
    paramsOut = [k1N', n1N', tN', k0N'];
    morphologyOut = mN';
    
end