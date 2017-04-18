% parameter extraction for clinical study embryos
% Livia Zarnescu Yanez 11-10-15
%

clear all;
warning off;
close all;
pNum = '033';
date = '9-18-16';
embryoNum = '8';

manualPip = 0;
manualCorner = 0;
manualMeasure = 0;
pressureUsed = 0.4;
filePath = 'C:\Users\Livia\Desktop\IVF\';

filePathRaw = [filePath 'Raw Data\Videos\Human\MECH' pNum];
if ~isempty(date)
    filePathRaw = [filePathRaw '-' date];
end

filePathProc = [filePath 'Processed Data\Human\MECH' pNum];
if ~exist(filePathProc, 'dir')
    mkdir(filePathProc);
end

fileNameSave = [filePathProc '\MECH' pNum '_E' embryoNum];
movpath = [filePathRaw '\MECH' pNum '-E' embryoNum '.avi'];
% movpath = [filePathRaw '\E' embryoNum '.avi'];

% some parameters ...
convFactor = 2.27; % pixels / micron on clinical system
cropVal = 1; % maximum num seconds to fit to model
secsToGet = 1.2;
startFrame = 28;
cannyThresh = .08;
pixelSag = 20;

% read in video
[newframes, frameRate] = ReadInVideo(movpath, secsToGet, startFrame, 1, 1);

% get pipette ROI
% make new pipette reference if one does not already exist
if manualPip
    [~, ROIframes] = getROI(newframes, manualCorner, cannyThresh);
%     figure, imshow(ROIframes(:,:,1));
%     [x,y] = ginput(2);
%     close all;
%     pipRefOpeningPixels = round(abs(y(2) - y(1)));
else
    ROIframes = GetPipetteROI(newframes, cannyThresh, NaN, filePathRaw, 0);
    load([filePathRaw '\pipRef.mat']);
end

t = 0:(1/frameRate):((size(ROIframes,3)-1)/frameRate);
clear newframes;

% make pipette size reference, convert opening size to microns
Fin = pressureUsed * 0.87 * 6895 * pi * ((pipRefOpeningPixels/(convFactor*2))*10^-6)^2; % pressure*area,

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
A = (aspiration_depth - pixelSag) * 10^-6 / convFactor; % convert from pixels to meters

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

save([fileNameSave '.mat'], 'xfine', 'yfit', ...
        'k0', 'k1', 'n0', 'n1', 'tau', 'F0', 'fval', 't', ...
        'aspiration_depth', 'A');



















