% load embryo video for measurement
% LZ 8-31-15

function [newframes, t, params] = LoadVideo(type, currDate, pipSize, ...
    filePathRaw, embryoNum, handles)

dateU = currDate;
dateUI = strfind(currDate, '-');
dateU(dateUI) = '_';
frameStartMult = .45;
cannyThresh = .35;

% get applied pressure from file
pressureFileName = [filePathRaw '\E' num2str(embryoNum) '_pressure.txt'];
f = fopen(pressureFileName);
pressureList = fscanf(f, '%f');
fclose(f);
pressureApplied = pressureList(end);

if isequal(type, 'Human')
    procFileName = ['aspiration_data_', dateU, '_human_E'];
    frameMult = 1.2;
    cropVal = 1;
    pipSizeRef = 229;
    Fin = pressureApplied * 6895 * pi * (35*pipSize/pipSizeRef*10^-6)^2; % pressure*area, 
    % 229 pixels is standard pipette opening
    adaptHist = 0;
elseif isequal(type, 'Mouse Oocyte')
    procFileName = ['aspiration_data_', dateU, '_E'];
    frameMult = .8;
    cropVal = .5;
    pipSizeRef = 126;
    Fin = pressureApplied * 6895 * pi * (20*pipSize/pipSizeRef*10^-6)^2; % pressure*area
    adaptHist = 1;
elseif isequal(type, 'Mouse Embryo')
    procFileName = ['aspiration_data_', dateU, '_E'];
    frameMult = .8;
    cropVal = .5;
    pipSizeRef = 126;
    Fin = pressureApplied * 6895 * pi * (20*pipSize/pipSizeRef*10^-6)^2; % pressure*area
    % 126 pixels is standard opening
    adaptHist = 1;
end

if strcmp(currDate,'11-19') || strcmp(currDate,'11-30')
    convFactor = 54;
    pipLarge = 0;
else
    convFactor = 128;
    pipLarge = 1;
end

% set path of video and time stamp files
movpath = [filePathRaw '\E' num2str(embryoNum) '.avi'];
timepath = [filePathRaw '\E' num2str(embryoNum) '.txt'];

if ~exist(movpath, 'file')
    errordlg('Error! Video file not found');
    return;
end

if ~exist(timepath, 'file')
    errordlg('Error! Text file (video timing) not found');
    return;   
end


secsToGet = -1;
startFrame = 1;
[newframes, frameRate] = ReadInVideo(movpath, secsToGet, startFrame, 0, 0);
t = 0:1/frameRate:((size(newframes,3)-1)/frameRate);

% % read in video and time stamp
% obj = VideoReader(movpath);
% fid = fopen(timepath);
% tcell = textscan(fid, '%f');
% fclose(fid);
% tlist = tcell{1,1};
% tlist = tlist(3:end)';
% frameRate = obj.FrameRate;
% 
% % read in frames starting just before aspiration
% currFirstFrame = 1;
% currLastFrame = min(obj.NumberOfFrames - 2, round(frameRate*(frameStartMult+2)));
% 
% 
% % read in some of the frames, define time vector
% frames = read(obj, [currFirstFrame currLastFrame]);
% tlist = tlist(currFirstFrame:currLastFrame);
% t = 0:1/frameRate:(currLastFrame-currFirstFrame)/frameRate;
% 
% % convert frames to grayscale and double format
% s = size(frames);
% numFrames = s(4);
% newframes = zeros(size(frames,1), size(frames,2), numFrames);
% for i = 1:numFrames
%     if s(3) > 1
%         % take each frame and make it grayscale
%         newframes(:,:,i) = double(rgb2gray(frames(:,:,:,i)))/255;
%     else
%         newframes(:,:,i) = double(frames(:,:,i))/255;
%     end
% end

imshow(newframes(:,:,round(frameStartMult*frameRate)), 'Parent', handles.MeasAxes);

% make params structure with random params to pass for later use
params.procFileName = procFileName;
params.frameMult = frameMult;
params.cropVal = cropVal;
params.pipSizeRef = pipSizeRef;
params.Fin = Fin;
params.adaptHist = adaptHist;
params.convFactor = convFactor;
params.pipLarge = pipLarge;
params.frameStartMult = frameStartMult;
params.cannyThresh = cannyThresh;
params.frameRate = frameRate;

