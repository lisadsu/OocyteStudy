% parameter extraction for human oocyte measurements
% based off depthDetectionAuto.m but more generalized
% no pressure file saved in clinical system, video only so relies on known
% pressure

% Livia Zarnescu Yanez
% 7-5-16

nameToSave = 'E1_GV'; 
manualPip = 0; % detect pipette edge based off its corner or not
manualCorner = 0; % manually select ROI
manualMeasure = 0;

addpath('C:\Users\Livia\Desktop\IVF\Code\EmbryoProject');
addpath('C:\Users\Livia\Desktop\IVF\Code\EmbryoProject\functions');
addpath('C:\Users\Livia\Desktop\IVF\Code\EmbryoProject\OocyteStudies');

basePath = 'C:/Users/Livia/Desktop/IVF/Raw Data/Videos/Human/HumanOocyteMeas';
cd(basePath);
[fileName, pathName] = uigetfile('*.avi');
pipSize = 30; % total pipette diameter
pressureApplied = .1;
frameStartMult = .45; % seconds at which to start video
frameMult = 1.2; % how many seconds of video to load
cannyThresh = .35;
convFactor = 2.27; % pixels / micron on clinical system
cropVal = 1;

% read in video
secsToGet = -1;
startFrame = 25;
[newframes, frameRate] = ReadInVideo([pathName fileName], secsToGet, startFrame, 1, 0);
currFrame = round(frameStartMult*frameRate);
newframes = newframes(:,:,1:round(frameMult*frameRate));

% get pipette ROI
% make new pipette reference if one does not already exist
if manualPip
    [~, ROIframes] = getROI(newframes, manualCorner, cannyThresh);
    
    if ~exist('pipRefOpeningPixels', 'var')
        figure(1);
        clf;
        imshow(ROIframes(:,:,1));
        h = helpdlg('Click on inner edges of pipette');
        pause(.5);
        close(h);
        [~,y] = ginput(2);
        pipRefOpeningPixels = abs(round(y(2) - y(1)));
    end
else
    ROIframes = GetPipetteROI(newframes, cannyThresh, NaN, pathName, 0);
    load([pathName 'pipRef.mat']);
end

clear newframes;
t = 0:(1/frameRate):((size(ROIframes,3)-1)/frameRate);
Fin = pressureApplied * .867 * 6895 * pi * ...
    ((pipRefOpeningPixels/(convFactor*2))*10^-6)^2; % pressure*area,



% meas offsetVal (initial depth offset inside pipette with only holding pressure)
figure(1);
clf;
imshow(ROIframes(:,:,1));
h = helpdlg('Click on initial cell depth');
pause(.5);
close(h);
[offsetVal,~] = ginput(1);


% measure parameters
if manualMeasure
       
    % depth detection
    % eventually automate this
    aspiration_depth = GetAspirationDepthManual(size(ROIframes,3),ROIframes);
    
    % adjust for varying vector lengths
    % chop off end if there is mismatch
    minLength = min(length(aspiration_depth), length(t));
    t = t(1:minLength);
    aspiration_depth = aspiration_depth(1:minLength);

else
    aspiration_depth = ...
        GetAspirationDepthAuto(size(ROIframes,3), ROIframes, 0);
end


% Fit to Model
% A is in METERS
% aspiration_depth is in PIXELS
close all;
t = t(t < cropVal);
aspiration_depth = aspiration_depth(1:length(t));
aspiration_depth = aspiration_depth(aspiration_depth > 0);
t = t(1:length(aspiration_depth));
A = (aspiration_depth - offsetVal) * 10^-6 / convFactor; % convert from pixels to meters

figure(5);
clf;

tauTryList = .02:.02:.2;
fValList = zeros(1,length(tauTryList));

for kk = 1:length(tauTryList)
    
    start_params(1) = .2; % k0
    start_params(2) = .2; % k1
    start_params(3) = tauTryList(kk);%tauStart; % tau
    start_params(4) = 5; % n1_inv (slope of creep)
    
    [xfine, yfit, k0, k1, n0, n1, F0, tau, fval] = KelvinFit3(t, A, Fin, 1, start_params);
    fValList(kk) = fval;
end

figure(5);
clf;
start_params(3) = tauTryList(fValList == min(fValList));
[xfine, yfit, k0, k1, n0, n1, F0, tau, fval] = KelvinFit3(t, A, Fin, 1, start_params);
fval
[k1 n1 tau k0]


save([basePath '/MeasuredData/' nameToSave '.mat'], 'xfine', 'yfit', ...
        'k0', 'k1', 'n0', 'n1', 'tau', 'F0', 'fval', 't', ...
        'aspiration_depth', 'A', 'offsetVal');












