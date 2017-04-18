%
% This function returns 3 parameters given an input structure with 
%   info about all the embryos. embryoNum is the element within the
%   structure to access. EmbryoInfo contains ROI, mask, x, y, r.
%
%   The parameters calculated are 1) the overall brightness of the
%   cytoplasm, 2) a profile of the brightness around the outline of each
%   cell, and 3) a profile of the average brightness going radially outward
%   from the center of each cell. 
%
function [brightParam, circParam, radParam] = getROIparams_A(embryoInfo, imageGroup, embryoNum)

groupNames = embryoInfo.groupNames;

cellROI = embryoInfo.(groupNames{imageGroup}).(['E' num2str(embryoNum)]).ROI;
cellMask = embryoInfo.(groupNames{imageGroup}).(['E' num2str(embryoNum)]).mask;
r = embryoInfo.(groupNames{imageGroup}).(['E' num2str(embryoNum)]).mR; 
xCenter = embryoInfo.(groupNames{imageGroup}).(['E' num2str(embryoNum)]).mX; % x coordinate of cell center
yCenter = embryoInfo.(groupNames{imageGroup}).(['E' num2str(embryoNum)]).mY; % y coordinate of cell center
figure(3);

% make mask with only cytoplasm pixels
cellCytoplasm = cellROI .* cellMask;

% extract intensity profile as a function of distance from center
% Use 20 bins. Initialize intensity profile.
rStepList = (0:.05:1) * r;
radParam = NaN*zeros(1, length(rStepList)-1);

% calculate current distance from center
[iC, jC] = meshgrid(1:size(cellROI,2), 1:size(cellROI,1));
distFromCenter = sqrt((iC - xCenter).^2 + (jC - yCenter).^2);

% calculate radial profile outward from center
for j = 1:length(rStepList)-1
    
    % curr radius range is from rStepList(j):rStepList(j+1)
    currDiskMask = (distFromCenter > rStepList(j)) & ...
        (distFromCenter < rStepList(j+1));
    intensityMask = cellROI .* currDiskMask .* cellMask;
    radParam(j) = mean(intensityMask(~isnan(intensityMask) & intensityMask > 0));
    
end

% calculate average brightness of cytoplasm pixels
brightParam = mean(cellCytoplasm(cellCytoplasm > 0));

% 
[xPath, yPath] = circlepoints(r*.9);
xPath = xPath + xCenter;
yPath = yPath + yCenter;
circParam = (improfile(cellROI, xPath, yPath, 100))';

% plot
switch embryoInfo.(groupNames{imageGroup}).(['E' num2str(embryoNum)]).group
    case 1
        plotColor = [1 .6 .6]; % Ab-soft: pale red
    case 2
        plotColor = [.8 .3 .3]; % Ab-mid: med red
    case 3
        plotColor = [.6 0 0]; % Ab-stiff: dark red
    case 4
        plotColor = [.6 .6 1]; % BAPTA-soft: pale blue
    case 5
        plotColor = [.3 .3 .8]; % BAPTA-mid: med blue
    case 6
        plotColor = [0 0 .6]; % BAPTA-stiff: dark blue
    case 7
        plotColor = [.6 1 .6]; % Control-soft: pale green
    case 8
        plotColor = [.3 .8 .3]; % Control-mid: mid green
    case 9
        plotColor = [0 .6 0]; % Control-stiff: dark green
    case 10
        plotColor = [1 1 .6]; % NonControl-soft: pale yellow
    case 11
        plotColor = [.8 .8 .3]; % NonControl-mid: mid yellow
    case 12
        plotColor = [.6 .6 0]; % NonControl-stiff: dark yellow
end


subplot(1,2,1);
hold on;
plot(radParam, 'color', plotColor);

subplot(1,2,2);
hold on;
p = polar(linspace(0, 2*pi, length(circParam)), circParam);
set(p, 'color', plotColor);


