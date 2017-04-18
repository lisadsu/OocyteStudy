% read in a video file and convert it to grayscale

function [framesOut, frameRate] = ReadInVideo(movpath, vidSecs, startFrame, flipside, flipvert)

if nargin < 4
    flipside = 0;
    flipvert = 0;
end

if nargin < 5
    flipvert = 0;
end

% read in video
obj = VideoReader(movpath);
frameRate = obj.FrameRate;

% get all or some
if vidSecs == -1
    totalFrames = floor(frameRate*obj.Duration);
else
    totalFrames = round(frameRate*vidSecs);
end

obj.CurrentTime = startFrame / obj.FrameRate;
lastFrame = round(min(obj.FrameRate*obj.Duration - 1, startFrame + totalFrames));
numFrames = lastFrame - startFrame + 1;


% convert to grayscale and double format
framesOut = zeros(obj.Width, obj.Height, numFrames);
%framesOut = zeros(obj.Height, obj.Width, numFrames);

for i = startFrame:lastFrame
    
    currFrame_pre = readFrame(obj);
    %flip
    currFrame = imrotate(currFrame_pre,-90);
    s = size(currFrame);
    
    if length(s) > 2
        % take each frame and make it grayscale
        if flipside && flipvert
            framesOut(:,:,i-startFrame+1) = rot90(double(rgb2gray(currFrame)),2)/255;
        elseif flipvert
            framesOut(:,:,i-startFrame+1) = flipud(double(rgb2gray(currFrame)))/255;
        elseif flipside
            framesOut(:,:,i-startFrame+1) = fliplr(double(rgb2gray(currFrame)))/255;
        else
            framesOut(:,:,i-startFrame+1) = double(rgb2gray(currFrame))/255;
        end
    else
        if flipside && flipvert
            framesOut(:,:,i-startFrame+1) = rot90(double(currFrame),2)/255;
        elseif flipvert
            framesOut(:,:,i-startFrame+1) = flipud(double(currFrame))/255;
        elseif flipside
            framesOut(:,:,i-startFrame+1) = fliplr(double(currFrame))/255;
        else
            framesOut(:,:,i-startFrame+1) = double(currFrame)/255;
        end
    end
    
%     framesOut(:,:,i-startFrame+1) = imrotate(framesOut(:,:,i-startFrame+1),-20, 'nearest', 'crop');
end





