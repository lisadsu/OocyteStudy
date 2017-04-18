% Go through all confocal images of cortical granules and extract params
% This code does not apply to microinjection experiments
% For re-analysis, just load the params
% Livia Zarnescu
% 1-27-15

% load struct with all 36 embryos, and their associated file paths
load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Embryo\cg_imaging\embryoInfoCG.mat');
nEmbryos = sum(embryoInfo.numEmbryosInImage);
newImDisplayIndex = [1 cumsum(embryoInfo.numEmbryosInImage) + 1];
embryoNumSum = [0 cumsum(embryoInfo.numEmbryosInImage)];

%%
    
% for each figure, select an ROI to analyze
for i = 1:nEmbryos
    
    % load and display image
    currFileName = embryoInfo.(['E' num2str(i)]).filePath{1};
    imageNum = find(i <= cumsum(numEmbryosInImage), 1, 'first');
%     embryoInfo.(['E' num2str(i)]).embryoNumInImage = i - embryoNumSum(imageNum);
    
    % if we have to display a new file
    if ~isempty(intersect(i, newImDisplayIndex))
        
        % calculate projection image
        imInfo = imfinfo(currFileName); % get dimensions of first image
        imageMatrix = zeros(imInfo(1).Height, imInfo(1).Width); % green channel
        
        for k = 1:length(imInfo)
            currImage = im2double(imread(currFileName, 'Index', k));
            imageMatrix = imageMatrix + currImage;
        end
        
        % take average and apply blur
        imageMatrix = imageMatrix / length(imInfo);
        imageMatrix = imfilter(imageMatrix, fspecial('gaussian', 10, 2));
%         imMatrixColor = zeros(size(imageMatrix, 1), size(imageMatrix, 2), 3);
%         imMatrixColor(:,:,2) = imageMatrix;
        
        % display image
        figure(1);
        clf;
        imagesc(imageMatrix);
        set(gca, 'fontsize', 14);
        
    end
    
 
    % select ROI, then calc "brightness" normalized by area
    figure(1);
%     title(['draw region around embryo ' num2str(i - embryoNumSum(imageNum))]);
%     coords = getrect;
%     coords = embryoInfo.(['E' num2str(i)]).coords;
%     ROI = imadjust(imageMatrix(round(coords(2):(coords(2)+coords(4))), ...
%         round(coords(1):(coords(1)+coords(3)))), [0 .05], [0 1]);
    ROI = embryoInfo.(['E' num2str(i)]).ROI;
    figure(2);
    clf;
    imshow(imadjust(ROI));
    
    % now calculate params!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % param for overall cytoplasm brightness, maybe brightness going
    % radially outward, and maybe brightness around outside of cell
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [brightParam, circParam, radParam] = getROIparams(embryoInfo, i);
    
    embryoInfo.(['E' num2str(i)]).brightParam = brightParam;
    embryoInfo.(['E' num2str(i)]).circParam = circParam;
    embryoInfo.(['E' num2str(i)]).radParam = radParam;
    embryoInfo.(['E' num2str(i)]).circStd = std(circParam);
    embryoInfo.(['E' num2str(i)]).circMean = mean(circParam);
    
    % get mask and save
    % ask user if they want to manually select mask or detect automatically
    % with circle hough transform
%     
%     maxRadius = mean(size(ROI))/2;
%     radiusRange = round(linspace(maxRadius*.8, maxRadius, 20))
%     [mask, mR, mX, mY] = getCellMask(ROI, radiusRange);
%     
%     embryoInfo.(['E' num2str(i)]).ROI = ROI;
%     embryoInfo.(['E' num2str(i)]).mask = mask;
%     embryoInfo.(['E' num2str(i)]).x = mX;
%     embryoInfo.(['E' num2str(i)]).y = mY;
%     embryoInfo.(['E' num2str(i)]).r = mR;
%     
    
    % save coords
%     embryoInfo.(['E' num2str(i)]).coords = coords;
    
    % brightness is average brightness of nonzero pixels squared
%     brightParam = squeeze(sum(sum(double((ROI(:,2,:)).*(ROI(:,2,:) > 0)) / sum(sum(ROI(:,2,:) > 0)))))
%     granuleParam = [granuleParam brightParam]
    
%     embryoInfo.(['E' num2str(i)]).brightParam = brightParam;
    
end

% embryoInfo.granuleParam = granuleParam;
% save('embryoInfoCG.mat', 'embryoInfo');


%% Make figures from params

circMeanList = [];
circStdList = [];
radMeanList = [];
brightParamList = [];
circDiffList = [];
% 
% go thru embryo struct, accumulate list of params to plot
% for i = 1:nEmbryos
% 
%     circMeanList = [circMeanList embryoInfo.(['E' num2str(i)]).circMean];
%     circStdList = [circStdList embryoInfo.(['E' num2str(i)]).circStd];
%     radMeanList = [radMeanList mean(embryoInfo.(['E' num2str(i)]).radParam)];
%     brightParamList = [brightParamList embryoInfo.(['E' num2str(i)]).brightParam]; 
%     circDiffList = [circDiffList (max(embryoInfo.(['E' num2str(i)]).circParam) ... 
%         - min(embryoInfo.(['E' num2str(i)]).circParam))/embryoInfo.(['E' num2str(i)]).circMean];
%     
% end

% embryoInfo.circMeanList = circMeanList;
% embryoInfo.circStdList = circStdList;
% embryoInfo.radMeanList = radMeanList;
% embryoInfo.brightParamList = brightParamList;
% embryoInfo.circDiffList = circDiffList;

% save('./data/embryoInfoCG.mat', 'embryoInfo');

circMeanList = embryoInfo.circMeanList;
circStdList = embryoInfo.circStdList;
circDiffList = embryoInfo.circDiffList;
radMeanList = embryoInfo.radMeanList;
brightParamList = embryoInfo.brightParamList;
viability = embryoInfo.viability;

colorVec = [0 0 .6; ... % nonviable
            0 .6 0];    % viable
        
figure, scatter3(circMeanList, brightParamList, circDiffList, ...
    100, colorVec(viability + 1,:), 'filled');
xlabel('x');
ylabel('y');
zlabel('z');

paramToPlot = 'circMeanList';

figure;
p1 = bar(.75, mean(eval([paramToPlot '(viability == 0)'])), .4, 'facecolor', colorVec(1,:));
hold on;
p1e = errorbar(.75, mean(eval([paramToPlot '(viability == 0)'])), ...
    std(eval([paramToPlot '(viability == 0)'])), 'color', 'k', 'linewidth', 2);

p2 = bar(1.25, mean(eval([paramToPlot '(viability == 1)'])), .4, 'facecolor', colorVec(2,:));
hold on;
p2e = errorbar(1.25, mean(eval([paramToPlot '(viability == 1)'])), ...
    std(eval([paramToPlot '(viability == 1)'])), 'color', 'k', 'linewidth', 2);

set(gca, 'fontsize', 14);
set(gca, 'xtick', [.75 1.25])
set(gca, 'xticklabel', {'Non-viable', 'Viable'});
ylabel('Parameter value (a.u.)');
title(sprintf('Average brightness \nnear cell membrane'));




%% Load struct with all embryo params and plot

% load struct with all 36 embryos, and their associated file paths
load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Embryo\cg_imaging\embryoInfoCG.mat');
nEmbryos = sum(embryoInfo.numEmbryosInImage);
newImDisplayIndex = [1 cumsum(embryoInfo.numEmbryosInImage) + 1];
    
for i = 1:nEmbryos
    
    % load and display image
    currFileName = embryoInfo.(['E' num2str(i)]).filePath{1};
    
    % if we have to display a new file
    if ~isempty(intersect(i, newImDisplayIndex))
        
        % calculate projection image
        imInfo = imfinfo(currFileName); % get dimensions of first image
        imageMatrix = zeros(imInfo(1).Height, imInfo(1).Width); % green channel
        
        for k = 1:length(imInfo)
            currImage = im2double(imread(currFileName, 'Index', k));
            imageMatrix = imageMatrix + currImage;
        end
        
        % take average and apply blur
        imageMatrix = imageMatrix / length(imInfo);
        imageMatrix = imfilter(imageMatrix, fspecial('gaussian', 10, 2));
        
        % display image
        %figure(1);
        %clf;
        %imagesc(imageMatrix);
        %set(gca, 'fontsize', 14);
        
    end
    
 
    % select ROI, then calc "brightness" normalized by area
    %figure(1);
    ROI = embryoInfo.(['E' num2str(i)]).ROI;
    %figure(2);
    %clf;
    %imshow(imadjust(ROI));
    
    % plot
    if embryoInfo.(['E' num2str(i)]).viability == 1
        plotColor = [0 .6 0];
    else
        plotColor = [0 0 .6];
    end
    
    figure(3);
    hold on;
    circParam = embryoInfo.(['E' num2str(i)]).circParam;
    p = plot(linspace(0, 2*pi, length(circParam)), circParam, 'linewidth', 2);
    set(p, 'color', plotColor);

    
end

set(gca, 'fontsize', 14);
xlim([0 2*pi]);
xlabel('Angle (rad)');
ylabel('Image Intensity (a.u.)');
h = get(gca, 'children');
legend([h(8) h(20)], 'viable', 'nonviable');
title('Image Intensity Around Cell Boundary');

% for (i = 1:length(h))
%     set(h(i), 'linewidth', 2);
% end



