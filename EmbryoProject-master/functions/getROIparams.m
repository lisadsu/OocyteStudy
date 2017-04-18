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
function [brightParam, circParam, radParam] = getROIparams(embryoInfo, embryoNum)

cellROI = embryoInfo.(['E' num2str(embryoNum)]).ROI;
cellMask = embryoInfo.(['E' num2str(embryoNum)]).mask;
r = embryoInfo.(['E' num2str(embryoNum)]).r; 
xCenter = embryoInfo.(['E' num2str(embryoNum)]).x; % x coordinate of cell center
yCenter = embryoInfo.(['E' num2str(embryoNum)]).y; % y coordinate of cell center
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
if embryoInfo.(['E' num2str(embryoNum)]).viability == 1
    plotColor = [0 0 1];
else
    plotColor = [1 0 0];
end

subplot(1,2,1);
hold on;
plot(radParam, 'color', plotColor);

subplot(1,2,2);
hold on;
p = polar(linspace(0, 2*pi, length(circParam)), circParam);
set(p, 'color', plotColor);


