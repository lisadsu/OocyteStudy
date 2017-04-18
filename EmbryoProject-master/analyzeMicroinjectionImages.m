% extract parameters from cortical granule images
% after microinjection experiments
% Livia Zarnescu
% 10-1-14

clear all;
close all;

%% 1. For each embryo, save tiff stack

% run this ONLY if TIFFs don't already exist
imageDirectory = 'C:\Users\Livia\Desktop\IVF\CG stain\microinjections\';
dirNames = {'Ab-soft\', 'Ab-mid\', 'Ab-stiff\', 'BAPTA-soft\CG-1\', ...
    'BAPTA-soft\CG-2\',  'BAPTA-mid\', 'BAPTA-stiff\', ...
    'Control-soft\', 'Control-mid\', 'Control-stiff\', ...
    'Non-control-soft\1\', 'Non-control-soft\2\', ...
    'Non-control-mid\', 'Non-control-stiff\'};
fileNames = {'soft', 'mid', 'stiff', 'soft-1', 'soft-2', 'mid', ...
    'stiff-BAPTA', 'soft', 'mid', 'stiff', 'noninjected-soft-1', ...
    'noninjected-soft-2', 'noninjected-mid', 'noninjected-stiff'};
tifNames = {'Ab-soft', 'Ab-mid', 'Ab-stiff', 'BAPTA-soft1', ...
    'BAPTA-soft2', 'BAPTA-mid', 'BAPTA-stiff', ...
    'Control-soft', 'Control-mid', 'Control-stiff', ...
    'Non-control-soft1', 'Non-control-soft2', 'Non-control-mid', ...
    'Non-control-stiff'};


for i = 1:length(dirNames)
    
    display(['Filename: ' tifNames(i)]);
    
    % get list of all tif images in directory and 
    % make empty matrix to hold them
    imageList = dir([imageDirectory dirNames{i}]);
    imageList = imageList(3:end);
    imageProps = imfinfo([imageDirectory dirNames{i} imageList(1).name]);
    imageMatrix = zeros(imageProps.Width, imageProps.Height, length(imageList));
    
    for j = 1:length(imageList)
        currFileBaseName = [imageDirectory dirNames{i} imageList(j).name];
        A = imread(currFileBaseName);
        A = im2double(A(:,:,2));
        imageMatrix(:,:,j) = A;
    end
    
    figure, imagesc(mean(imageMatrix,3));
    colorbar;
    caxis([0 .15]);
    pause;
    
    WriteTiffStack(imageMatrix, [imageDirectory 'TIFS\' tifNames{i} '.tif']);
    
end


%% 2. Make struct containing ROIs for each embryo

% build struct
% sub-structs for group, and embryo within each group
fileToSave = 'C:/Users/Livia/Desktop/IVF/Processed Data/Mouse Embryo/microinjection_data/embryoInfoMicroinjection.mat';
groupNames = {'Ab_soft', 'Ab_mid', 'Ab_stiff', ...
    'BAPTA_soft', 'BAPTA_mid', 'BAPTA_stiff', ...
    'Control_soft', 'Control_mid', 'Control_stiff', ...
    'NonControl_soft', 'NonControl_mid', 'NonControl_stiff'};
embryosInImage = [4 7 5 1 1 8 5 5 5 6 1 1 6 7];
embryosInGroup = [4 7 5 2 8 5 5 5 6 2 6 7];
imageGroup = [1 2 3 4 4 5 6 7 8 9 10 10 11 12];
eNumOffset = [0 0 0 0 1 0 0 0 0 0 0 1 0 0]; % for groups with multiple TIFS
embryoInfo = struct;

for i = 1:length(groupNames)
    embryoInfo.(groupNames{i}) = struct;
    for j = 1:sum(embryosInImage(imageGroup == i))
        embryoInfo.(groupNames{i}).(['E' num2str(j)]) = struct;
        
        if i == 1 && j == 1
            embryoInfo.(groupNames{i}).(['E' num2str(j)]).tifName = tifNames{1};
        else
            embryoInfo.(groupNames{i}).(['E' num2str(j)]).tifName = tifNames{i+1};
        end
    end
end

embryoInfo.embryosInImage = embryosInImage;
embryoInfo.groupNames = groupNames;
embryoInfo.imageGroup = imageGroup;
embryoInfo.tifNames = tifNames;
embryoInfo.imageDirectory = imageDirectory;
embryoInfo.embryosInGroup = embryosInGroup;
embryoInfo.eNumOffset = eNumOffset;

save(fileToSave, 'embryoInfo');

%% Read in Tiffs and extract ROIs

load(fileToSave);
embryosInImage = embryoInfo.embryosInImage;
groupNames = embryoInfo.groupNames;
eNumOffset = embryoInfo.eNumOffset;

for i = 1:length(embryoInfo.tifNames)
    
    currFileName = [embryoInfo.imageDirectory '\TIFS\' embryoInfo.tifNames{i} '.tif'];
    
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
    figure(1);
    clf;
    imagesc(imageMatrix);
    set(gca, 'fontsize', 14);
    
    imageGroup = embryoInfo.imageGroup(i);
    
    % select ROI and save
    for j = 1:embryosInImage(i)
        
        figure(1);
        title(['draw region around embryo ' num2str(j)]);
        
        % second image will have second embryo from first group
        eNum = j + eNumOffset(i);
        display(['Tiff # ' num2str(i), ', embryo #' num2str(eNum)]);
            
        coords = getrect;
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).coords = coords;
        ROI = imageMatrix(round(coords(2):(coords(2)+coords(4))), ...
            round(coords(1):(coords(1)+coords(3))));
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).ROI = ROI;
        
        figure(2);
        clf;
        imshow(imadjust(ROI));
        
    end
end

save(fileToSave, 'embryoInfo');

%% Detect cell boundary and cytoplasm

close all;
load(fileToSave);
embryosInImage = embryoInfo.embryosInImage;
groupNames = embryoInfo.groupNames;
eNumOffset = embryoInfo.eNumOffset;

for i = 1:length(embryoInfo.tifNames)
    
    imageGroup = embryoInfo.imageGroup(i);
    
    % select ROI and save
    for j = 1:embryosInImage(i)

        figure(1);
        title(['draw region around embryo ' num2str(j)]);
        
        % second image will have second embryo from first group
        eNum = j + eNumOffset(i);
        display(['Tiff # ' num2str(i), ', Group # ' ...
            num2str(imageGroup) ', embryo #' num2str(eNum)]);
            
        ROI = embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).ROI;
        
        figure(2);
        clf;
        imshow(imadjust(ROI));
        
        maxRadius = mean(size(ROI))/2;
        radiusRange = round(linspace(maxRadius*.7, maxRadius, 20));
        [mask, mR, mX, mY] = getCellMask(ROI, radiusRange);    
        
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).mask = mask;
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).mR = mR;
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).mX = mX;
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).mY = mY;
                
    end
end

save(fileToSave, 'embryoInfo');


%% Extract cell params

close all;
load(fileToSave);
embryosInGroup = embryoInfo.embryosInGroup;
groupNames = embryoInfo.groupNames;

for i = 9%1:length(embryosInGroup)
    
    % select ROI and save
    for j = 1:embryosInGroup(i)
        
        ROI = embryoInfo.(groupNames{i}).(['E' num2str(j)]).ROI;
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).group = i;
        
        figure(1);
        clf;
        imshow(imadjust(ROI));
        
         % get mean brightness of representative ROI
        coords = getrect;
        subROI = ROI(round(coords(2):(coords(2)+coords(4))), ...
            round(coords(1):(coords(1)+coords(3))));
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).subROI = subROI;
        
        % now calculate params!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % param for overall cytoplasm brightness, maybe brightness going
        % radially outward, and maybe brightness around outside of cell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [brightParam, circParam, radParam] = getROIparams_A(embryoInfo, i, j);

        embryoInfo.(groupNames{i}).(['E' num2str(j)]).brightParam = brightParam;
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).circParam = circParam;
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).radParam = radParam;
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).circStd = std(circParam);
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).circMean = mean(circParam);
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).subROImean = mean(mean(subROI));
        
    end
end

save(fileToSave, 'embryoInfo');


%% Make figures from params

load(fileToSave);
embryosInGroup = embryoInfo.embryosInGroup;
groupNames = embryoInfo.groupNames;
circMeanList = cell(1,12);
circStdList = cell(1,12);
radMeanList = cell(1,12);
brightParamList = cell(1,12);
circDiffList = cell(1,12);
ROImeanList = cell(1,12);
groupList = cell(1,12);
figure(1);
clf;
hold on;

for i = 1:length(embryosInGroup)
    
    % accumulate param vectors
    for j = 1:embryosInGroup(i)
        
        groupList{i} = [groupList{i} embryoInfo.(groupNames{i}).(['E' num2str(j)]).group];
        circMeanList{i} = [circMeanList{i} embryoInfo.(groupNames{i}).(['E' num2str(j)]).circMean];
        circStdList{i} = [circStdList{i} embryoInfo.(groupNames{i}).(['E' num2str(j)]).circStd];
        radMeanList{i} = [radMeanList{i} mean(embryoInfo.(groupNames{i}).(['E' num2str(j)]).radParam)];
        brightParamList{i} = [brightParamList{i} embryoInfo.(groupNames{i}).(['E' num2str(j)]).brightParam];
        circDiffList{i} = [circDiffList{i} (max(embryoInfo.(groupNames{i}).(['E' num2str(j)]).circParam) ...
            - min(embryoInfo.(groupNames{i}).(['E' num2str(j)]).circParam))/...
            embryoInfo.(groupNames{i}).(['E' num2str(j)]).circMean];
        ROImeanList{i} = [ROImeanList{i} embryoInfo.(groupNames{i}).(['E' num2str(j)]).subROImean];
        
        if mod(i, 3) == 2
            currColor = [0 .6 0];
            h1 = plot(linspace(0, 2*pi, 100), ...
                smooth(embryoInfo.(groupNames{i}).(['E' num2str(j)]).circParam), ...
                'color', currColor, 'linewidth', 2);
        else 
            currColor = [0 0 .6];
            h2 = plot(linspace(0, 2*pi, 100), ...
                smooth(embryoInfo.(groupNames{i}).(['E' num2str(j)]).circParam), ...
                'color', currColor, 'linewidth', 2);
        end
        
        
    
    end
end

set(gca, 'fontsize', 14);
xlabel('angle (rad)');
ylabel('image intensity (a.u.)');
title('Image Intensity Around Cell Boundary');
xlim([0 2*pi]);
legend([h1 h2], 'viable', 'nonviable');

embryoInfo.groupList = groupList;
embryoInfo.circMeanList = circMeanList;
embryoInfo.circStdList = circStdList;
embryoInfo.radMeanList = radMeanList;
embryoInfo.brightParamList = brightParamList;
embryoInfo.circDiffList = circDiffList;
embryoInfo.ROImeanList = ROImeanList;
save(fileToSave, 'embryoInfo');









