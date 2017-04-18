function auto_depth = GetAspirationDepthAuto(framesToGet, ...
    ROIframes, adaptHist, extraFig)

% Get user to indicate first frame where cell appears
% click on cell edge for first 5 frames

frameNum = 1;
currCellFrame = 1;
foundCell = 0;
manual_depth = zeros(1,5);

if nargin < 4 || ~ishandle(extraFig)
    fig = figure;
else
    fig = extraFig;
end

while ~foundCell || currCellFrame < 6
    
    % show current frame
    figure(fig);
    imshow(ROIframes(:,:,frameNum), 'InitialMagnification', 600);
    
    % make sure frame is in same position on screen each time
    if frameNum == 1
        figPos = get(fig, 'Position');
    else 
        set(fig, 'Position', figPos);
    end
    
    % ask user to click on edge of cell
    % xb will be <0 until cell appears
    [xb yb] = ginput(1);
    
    % if value is greater than 0, record cell edge for first 5 frames
    if xb > 0
        manual_depth(currCellFrame) = xb;
        if ~foundCell
            foundCell = 1;
        end
        currCellFrame = currCellFrame + 1;
    end
    
    frameNum = frameNum+1;
    
end

frameNum = frameNum - 1;

% now frameNum will have the value of the first frame to start taking
% automated measurements

%% Get a rectangular region around last manual depth position in center of pipette

% grab box centered around yb, xb
% center of pipette in y and most recent xb in x

rSum = smooth(sum(ROIframes(:,:,1),2),20);
% figure, plot(rSum);
rRange = max(rSum) - min(rSum);
pipThresh = graythresh((rSum - min(rSum))/rRange)*rRange + min(rSum);

currInd = 1;
startInd = 1;
endInd = 1;
startP = 0;
endP = 0;

while ~endP
    
    % watch for when value goes above threshold
    if rSum(currInd) > pipThresh && ~startP
        startP = 1;
        startInd = currInd;
    elseif rSum(currInd) < pipThresh && startP
        endP = 1;
        endInd = currInd;
    end
        
    currInd = currInd + 1;
    
end
    
pipMidY = mean([startInd endInd]);
pipMidX = xb;

currY = yb;%pipMidY
currX = pipMidX;

scaleFactor = 5;
yRange = 10;
xRangeM = 10;
xRangeP = 10;
RF = imresize(medfilt2(imadjust(abs(ROIframes(:,:,frameNum)-ROIframes(:,:,1))),[3 3]), scaleFactor);
RF(RF < 0) = 0;
RF(RF > 1) = 1;

% get box from previous frame
firstBox = RF(round(scaleFactor*(currY-yRange):scaleFactor*(currY+yRange)),...
    round(scaleFactor*(currX-xRangeM):scaleFactor*(currX+xRangeP)));
firstBox = imadjust(firstBox);

%% Now match box to next frame, and repeat

auto_depth = zeros(1, framesToGet-5);
figure(fig);
clf;
currBox = firstBox;
% mN = min(min(medfilt2(abs(ROIframes(:,:,frameNum+1)-ROIframes(:,:,1)),[3 3])));
% mX = max(max(medfilt2(abs(ROIframes(:,:,frameNum+1)-ROIframes(:,:,1)),[3 3])));
% 
% fileName = 'C:\Users\Livia\Desktop\IVF\Processed Data\Human analysis\11-8-13 analysis human\aspirationVid.avi';
% RFcolor = zeros(size(ROIframes,1)*scaleFactor, size(ROIframes,2)*scaleFactor, 3);
% writerObj = VideoWriter(fileName);
% open(writerObj);

if (framesToGet + frameNum) > size(ROIframes,3)
    maxI = framesToGet - frameNum - 1;
else 
    maxI = framesToGet;
end

for i = 1:maxI
    
    % match currBox to next pipette frame
    if adaptHist
        RF = imresize(adapthisteq(medfilt2(abs(ROIframes(:,:,frameNum+i) - ...
            ROIframes(:,:,1)),[3 3])), scaleFactor);
    else
        RF = imresize(imadjust(medfilt2(abs(ROIframes(:,:,frameNum+i) - ...
            ROIframes(:,:,1)),[3 3])), scaleFactor);
    end
    RF(RF < 0) = 0;
    RF(RF > 1) = 1;
    
    [~,I_NCC,~] = template_matching(currBox, ...
        RF(round(scaleFactor*(currY-yRange):scaleFactor*(currY+yRange)),:));
    I_NCC = I_NCC(round(size(I_NCC,1)/2-scaleFactor*2):round(size(I_NCC,1)/2+scaleFactor*2),:);
    I_NCC = abs(I_NCC - .5);
    [yBest,xBest] = find(I_NCC == max(I_NCC(:)));
    xBest = xBest/scaleFactor - (xRangeM + xRangeP)/2 + xRangeM;
    yBest = yBest/scaleFactor;
    
    figure(fig);
    pause(.1);
    imshow(RF, 'InitialMagnification', 400);
    hold on;
    line(scaleFactor*[xBest xBest], scaleFactor*[currY-80 currY+80], 'color', [1 1 0], 'linewidth', 2);
    hold off;
    
%     % save tiff file with current frame
%     RFcolor(:,:,1) = RF;
%     RFcolor(:,:,2) = RF;
%     RFcolor(:,:,3) = RF;
%     RFcolor(round(scaleFactor*(currY-80):scaleFactor*(currY+80)),...
%         round(scaleFactor*xBest-1:scaleFactor*xBest+1),1) = 1;
%     RFcolor(round(scaleFactor*(currY-80):scaleFactor*(currY+80)),...
%         round(scaleFactor*xBest-1:scaleFactor*xBest+1),2) = 1;
%     RFcolor(round(scaleFactor*(currY-80):scaleFactor*(currY+80)),...
%         round(scaleFactor*xBest-1:scaleFactor*xBest+1),3) = 0;
% 
%     writeVideo(writerObj, RFcolor);
%     
    %     p1 = [currY-80 xBest-50]
    %     p2 = [currY xBest]
    %     p3 = [currY+80 xBest-50]
    %     plotCircle(p1, p2, p3);
    

    
    
    hold off;
    
    auto_depth(i) = xBest;
    
    currBox = RF(round(scaleFactor*(currY-yRange):scaleFactor*(currY+yRange)),...
        round(scaleFactor*(xBest-xRangeM):scaleFactor*(xBest+xRangeP)));
    
    currBox = imadjust(currBox);
    
end

auto_depth = [manual_depth auto_depth];
% figure, plot(auto_depth);




