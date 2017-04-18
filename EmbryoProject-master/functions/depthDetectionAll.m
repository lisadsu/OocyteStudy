% Livia Zarnescu
% 8-14-13
% general file for detecting aspiration depth from videos

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FILL IN TYPE AND DATE HERE %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type = 'mouse oocyte'; % other values = "human", "mouse embryo" or "mouse oocyte"
currDate = '1-23-14';
dateU = currDate;
dateUI = strfind(currDate, '-');
dateU(dateUI) = '_';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% MODIFY THESE AS NECESSARY %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

embryoNum = 1;
manualCorner = 0;
manualMeasure = 0;
frameStartMult = .95;
cannyThresh = .35;

filePath = 'C:\Users\Livia\Desktop\IVF\';

if isequal(type, 'human')
    filePath2 = ['Raw Data\Videos\Human\human tests ', currDate, '\'];
    filePath3 = ['Processed Data\Human analysis\', currDate, ...
        ' analysis human\'];
    filePath4 = ['aspiration_data_', dateU, '_human_E'];
    frameMult = 1.2;
    cropVal = 1;
    tauStart = .06;
    Fin = .347 * 6895 * pi * (35*10^-6)^2; % pressure*area
elseif isequal(type, 'mouse oocyte')
    filePath2 = ['Raw Data\Videos\Mouse Oocytes\videos ', currDate, '\'];
    filePath3 = ['Processed Data\Mouse oocyte analysis\', currDate, ...
        ' analysis\'];
    filePath4 = ['aspiration_data_', dateU, '_E'];
    frameMult = .8;
    cropVal = .5;
    tauStart = .03;
    Fin = .167 * 6895 * pi * (20*10^-6)^2; % pressure*area
elseif isequal(type, 'mouse embryo')
    filePath2 = ['Raw Data\Videos\Mouse Embryos\videos ', currDate, '\'];
    filePath3 = ['Processed Data\Mouse embryo analysis\', currDate, ...
        ' analysis\'];
    filePath4 = ['aspiration_data_', dateU, '_E'];
    frameMult = .8;
    cropVal = .5;
    tauStart = .03;
    Fin = .347 * 6895 * pi * (20*10^-6)^2; % pressure*area
end

if strcmp(currDate,'11-19') || strcmp(currDate,'11-30')
    convFactor = 54;
    pipLarge = 0;
else
    convFactor = 108;
    pipLarge = 1;
end

% set path of video and time stamp files
movpath = [filePath filePath2 'E' num2str(embryoNum) '.avi'];
timepath = [filePath filePath2 'E' num2str(embryoNum) '.txt'];

% read in video and time stamp
obj = mmreader(movpath);
fid = fopen(timepath);
tcell = textscan(fid, '%f');
fclose(fid);
tlist = tcell{1,1};
tlist = tlist(3:end)';
frameRate = obj.FrameRate;

% read in frames starting just before aspiration
% make detection automatic in the future maybe?
currFirstFrame = 1;
currLastFrame = round(frameRate*3);
totalFrames = round(frameRate*frameMult);

% read in some of the frames, define time vector
frames = read(obj, [currFirstFrame currLastFrame]);
tlist = tlist(currFirstFrame:currLastFrame);
% t = cumsum(tlist);
% t = t - min(t);
% t = t/1000;
t = 0:1/frameRate:(currLastFrame-currFirstFrame)/frameRate;

% convert frames to grayscale and double format
s = size(frames);
numFrames = s(4);
newframes = zeros(size(frames,1), size(frames,2), numFrames);
for i = 1:numFrames
    if length(s) == 4
        % take each frame and make it grayscale
        newframes(:,:,i) = double(rgb2gray(frames(:,:,:,i)))/255;
    else
        newframes(:,:,i) = double(frames(:,:,i))/255;
    end
end

clear frames;
framesToGet = round(frameRate*frameMult);

% choose ROI starting at pipette opening
% bottom bound should not be below bottom of pipette
[coords, ROIframes] = getROI(newframes, manualCorner, cannyThresh);
sROI = size(ROIframes);
clear newframes;
        
if manualMeasure
    
    numAvgs = 1;
    
    % average manual measurements numAvgs times
    aspiration_vals = zeros(1, length(tlist));
    
    for ii = 1:numAvgs

        % depth detection
        % eventually automate this
        manual_depth = GetAspirationDepthManual(framesToGet, ...
            ROIframes(:,:,round(frameStartMult*frameRate):end));
        
        % adjust for varying vector lengths
        % chop off end if there is mismatch
        minLength = min(size(manual_depth,2), size(t,2));
        t = t(1:minLength);
        minLength = min(size(aspiration_vals,2),size(manual_depth,2));
        aspiration_vals = aspiration_vals(1:minLength) + manual_depth(1:minLength);
                
        ii
        pause(2);
        
    end
    
    aspiration_depth = aspiration_vals / numAvgs;

else
    aspiration_depth = ...
        GetAspirationDepthAuto(framesToGet, ...
        ROIframes(:,:,round(frameStartMult*frameRate):end));
end


% Fit to Model
close all;
t = t(t < cropVal);
aspiration_depth = aspiration_depth(1:length(t));
A = aspiration_depth * 40 * 10^-6 / convFactor; % convert from pixels to meters

% start_params(1) = .2; % k0
% start_params(2) = .2; % k1
% start_params(3) = .03;%tauStart; % tau
% start_params(4) = 5; % n1_inv (slope of creep)
% 
% [xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 1, start_params);

figure(5);
clf;

tauTryList = .02:.02:.18;
fValList = [];

for kk = 1:length(tauTryList)
    
    start_params(1) = .2; % k0
    start_params(2) = .2; % k1
    start_params(3) = tauTryList(kk);%tauStart; % tau
    start_params(4) = 5; % n1_inv (slope of creep)
    
    [xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 1, start_params);
    fValList = [fValList fval];
end

figure(5);
clf;
start_params(3) = tauTryList(fValList == min(fValList));
[xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 1, start_params);
fval
[k1 n1 tau k0]
% ylim([1*10^-5 max(yfit)*1.2]);


% make sure directories to save data exist, and then save
% if auto measure, compare to previous manual measurement
if manualMeasure
    
    % save params if good
    if ~exist([filePath filePath3], 'dir')
        mkdir([filePath filePath3]);
    end
    
    save([filePath filePath3 filePath4 num2str(embryoNum) ...
        '.mat'], 'xfine', 'yfit', ...
        'k0', 'k1', 'n0', 'n1', 'tau', 'F0', 'fval', 't', ...
        'aspiration_depth', 'A');

else
    
    % save params if good
    if ~exist([filePath filePath3 'AutoMeasure\'], 'dir')
        mkdir([filePath filePath3], 'AutoMeasure');
    end
    
    save([filePath filePath3 'AutoMeasure\' filePath4 num2str(embryoNum) ...
        '.mat'], 'xfine', 'yfit', ...
        'k0', 'k1', 'n0', 'n1', 'tau', 'F0', 'fval', 't', ...
        'aspiration_depth', 'A');
    
    pause;
    % compare automatic measurement to previously done manual one to see if
    % it's accurate
    load([filePath filePath3 filePath4 num2str(embryoNum) '.mat']);
    [xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 1, start_params);
    
end



%% keep only relevant variables or replot

filePath = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse embryo analysis\8-9-13 analysis\aspiration_data_8_9_13_E';
n = [4 11 12 15 16 17 18 20 22 23 28 29 31 35 38 39 43 45 47 53 32 41 46 52 56];

for iNum = 1:length(n)
    
    iNum
    if exist([filePath num2str(iNum) '.mat']) > 0
        
        load([filePath num2str(iNum) '.mat']);
        
        [k1 n1 tau k0]

        aspiration_depth = aspiration_depth - 5;
        A = aspiration_depth * 40 * 10^-6 / convFactor; % convert from pixels to meters
        
        figure(5);
        clf;
        
        tauTryList = .02:.02:.18;
        fValList = [];
        
        for kk = 1:length(tauTryList)
            
            start_params(1) = .2; % k0
            start_params(2) = .2; % k1
            start_params(3) = tauTryList(kk);%tauStart; % tau
            start_params(4) = 5; % n1_inv (slope of creep)
            
            [xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 1, start_params);
            fValList = [fValList fval];
        end
        
        figure(5);
        clf;
        start_params(3) = tauTryList(fValList == min(fValList));
        [xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 1, start_params);
        fval
        [k1 n1 tau k0]
        % ylim([1*10^-5 max(yfit)*1.2]);


        save([filePath ...
            num2str(iNum) '.mat'], 'xfine', 'yfit', ...
            'k0', 'k1', 'n0', 'n1', 'tau', 'F0', 'fval', 't', ...
            'aspiration_depth', 'A');
    end
    
end



