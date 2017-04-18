%
%   this function outputs aspiration depth vector as selected manually
%   by user
%
%   numFrames = max number of frames in which to detect aspiration depth
%   ROIframes = video of aspiration cropped to start at pipette opening
%
%   manual_depth = aspiration depth in pixels
%

function manual_depth = GetAspirationDepthManual(numFrames, ROIframes, extraFig)

frameNum = 1;
firstRun = 1;
manual_depth1 = [0];

if nargin < 3 || ~ishandle(extraFig)
    fig = figure;
else
    fig = extraFig;
end

while frameNum < numFrames
    
    if firstRun
        frameNum = frameNum + 1;
        imshow(ROIframes(:,:,frameNum), 'InitialMagnification', 600);
        %truesize(fig, 5*size(ROIframes(:,:,1)));
        figPos = get(fig, 'Position');
    else
        frameNum = frameNum + 1;
        imshow(ROIframes(:,:,frameNum), 'InitialMagnification', 600);
        %truesize(fig, 5*size(ROIframes(:,:,1)));
        set(fig, 'Position', figPos);
    end
    [xb yb] = ginput(1);
    manual_depth1 = [manual_depth1 xb];
    firstRun = 0;
    
    
end

first_pos = find(manual_depth1 > 2, 1);
manual_depth = manual_depth1(first_pos:end);