% make ROI template


function [] = MakeRoiTemplate(filePathRaw, cannyThresh, displayFig, singleFrame)

if nargin < 3 || ~ishandle(displayFig)
    displayFig = figure;
end

if nargin < 4
    vidFileRef = uigetfile([filePathRaw '\*.avi'], 'Select video file');
    obj = VideoReader([filePathRaw '\' vidFileRef]);
    singleFrame = read(obj, 1);

    if length(size(singleFrame)) > 2
        singleFrame = rgb2gray(singleFrame);
    end
end

BW = edge(singleFrame, 'canny', cannyThresh);
figure(displayFig);
clf;
imshow(BW);
title('draw region with top corner of pipette');
coord = getrect;
coords =[coord(1) coord(1)+coord(3) coord(2) coord(2)+coord(4)];

% given rectangle with top corner of pipette, extract coordinates of edge
% line and return "corner" point
[xout yout] = returnCorner(BW, coords);
figure(displayFig);
clf;
imshow(BW);
truesize(displayFig, 2*size(singleFrame(:,:)));
title('choose bottom right bound of ROI');
[xb yb] = ginput(1);
ROIframe = singleFrame(round(yout-20:yb+20), round(xout-10:xb));

figure(displayFig);
clf;
imshow(ROIframe);

h = helpdlg('Click on inner edges of pipette');
pause(.5);
close(h);
[x,y] = ginput(2);
pipRefOpeningPixels = abs(round(y(2) - y(1)));

% save .mat file (ROI is offset by 10 pixels in x)
save([filePathRaw '\pipRef.mat'], 'ROIframe', 'pipRefOpeningPixels');

