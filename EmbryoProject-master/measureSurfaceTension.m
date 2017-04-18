% measure surface tension of oocytes
% LZ 10-2-15

% E8,E9.E10 10-8-15 is great example of stress/strain characteristic curve 
date = '10-8-15';
fileDir = ['C:\Users\Livia\Desktop\IVF\Raw Data\Videos\Mouse Oocyte\videos ' date];
procDir = ['C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\' date ' analysis'];
type = 'Zona'; 
eNum = 10;
embryoNum = 'E10';
convToMicrons = 40/128;
% colorList = [0 0 .6; 0 .6 0; .6 0 0; 0 .6 .6; .6 0 .6];

if exist([procDir '\' type '\all' type '.mat'], 'file')
    load([procDir '\' type '\all' type '.mat']);
end

pressureList = [];
depthList = [];

% fixed scale conversion (screenshot to actual pixels)
pipCurrOpeningPixels = 64;
pipRefOpeningPixels = 132;
imScale = pipCurrOpeningPixels / pipRefOpeningPixels;
cellSizeMeasured = 0;

while (true)
    
    [fileName, pathName] = uigetfile([fileDir '\' embryoNum '*.jpg']);
    
    if fileName == 0
        break;
    end
    
    display(['Now measuring ' fileName]);
    
    A = rgb2gray(imread([pathName fileName]));
    f1 = figure(1);
    clf;
    imshow(A);
    movegui(f1, [-100, -100]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % record current pressure in psi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pressure = input('Enter pressure (psi): ');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % detect pipette edge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check if pipRef exists in fileDir
    if ~exist([fileDir '\pipRef.mat'], 'file')
        % save template of original
        MakeRoiTemplate(fileDir, .35);
    end
    
    % template matching from mat file
    load([fileDir '\pipRef.mat']);

    % find difference in scale (optional)
    figure(f1);
    h = helpdlg('Select ROI pipette only');
    pause(.5); close(h);
    coord = getrect;
    coords = round([coord(1) coord(1)+coord(3) coord(2) coord(2)+coord(4)]);
    Apip = A(coords(3):coords(4), coords(1):coords(2));
    Acrop = A((coords(3)-20):(coords(4)+20), (coords(1)-20):(coords(2)+20));
    ROIframeScaled = imresize(ROIframe, imScale);

    % now do template matching on scaled version 
    [~,I_NCC,~] = template_matching(ROIframeScaled, Acrop);
    [yBest,xBest] = find(I_NCC == max(I_NCC(:)));
    yCrop = round(yBest - size(ROIframeScaled,1)/2);
    xCrop = round(xBest - size(ROIframeScaled,2)/2);
    AcropPip = Acrop((yCrop+10):(yCrop + size(ROIframeScaled,1)-10), (xCrop+10):end, :);
    
    figure(f1);
    clf;
    imshow(AcropPip, 'InitialMagnification', 500);
    movegui(f1, [-100,-100]);
    
    % finally, select the depth into the pipette for the current pressure
    h = helpdlg('Select aspiration depth into pipette');
    pause(.5); close(h);
    [x, y] = ginput(1);
    
    if x < 0
        x = 0;
    end
    
    if ~cellSizeMeasured
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % also measure cell size
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f1 = figure(1);
        clf;
        imshow(A);
        movegui(f1, [-100, -100]);
        h = helpdlg('Measure cell size twice');
        pause(.5); close(h);
        diamList = zeros(1,2);
        
        for j = 1:2
            [xx,yy] = ginput(2);
            d = sqrt((xx(1) - xx(2))^2 + (yy(1) - yy(2))^2);
            diamList(j) = d/imScale*convToMicrons;
        end
        
        cellCurr = mean(diamList);
        
        cellSizeMeasured = 1;
    end
    
    % save pressure and depth (convert back to pixels first)
    pressureList = [pressureList pressure];
    depthList = [depthList x/imScale*convToMicrons];
    
end

depthAll{eNum} = depthList;
pressureAll{eNum} = pressureList;
cellSize(eNum) = cellCurr;

% plot
% close all;
figure(5);
clf;
scatter(pressureList - min(pressureList), depthList - min(depthList), ...
    150, 'marker', '+', 'linewidth', 3, ...
    'markeredgecolor', [.3 .7 .7]);
set(gca, 'fontsize', 14);
xlabel('pressure (psi)');
ylabel('depth (\mum)');
title('Embryo mechanical properties');
xlim([-.02 max(pressureList - min(pressureList)) + .02]);
ylim([-1 max(depthList - min(depthList)) + 1]);

% % polyfit
% p = polyfit(pressureList, depthList, 1);
% depthFit = polyval(p, pressureList);
% 
% hold on;
% plot(pressureList, depthFit, 'linewidth', 2, 'color', 'k');

% depthList is in MICRONS
if ~isempty(depthList) 
    save([procDir '\' type '\' embryoNum '_' lower(type) '.mat'], ...
        'pressureList', 'depthList', 'cellCurr');
    if exist('cellSize', 'var')
        save([procDir '\' type '\all' type '.mat'], 'depthAll', 'pressureAll', 'cellSize');
    else
        save([procDir '\' type '\all' type '.mat'], 'depthAll', 'pressureAll');
    end
end



%% Measure surface tension

% load settings
type = 'Membrane';
embryoNums = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', ...
    '21', '22', '23', '24', '25', '26', '27', '28', '29', '30'};

if exist([procDir '\' type '\all' type '.mat'])
    load([procDir '\' type '\all' type '.mat']);
end

f6 = figure(6);
clf;
set(gca, 'fontsize', 14);
xlabel('stress (psi)');
ylabel('aspiration depth (\mum)');
title('membrane measurement');
% colorList = [0 0 .6; 0 .6 0; .6 0 0; 0 .6 .6; .6 0 .6];

% cellSize = zeros(1, length(pressureAll));
% surfTension = zeros(1, length(pressureAll));
% eqPressure = zeros(1, length(pressureAll));
plotHandles = [];

for i = 1:length(embryoNums)
    
    if ~isempty(pressureAll{i})
        
        hold on;
        subplot(4,5,i);
        h = scatter(pressureAll{i} - min(pressureAll{i}), ...
            depthAll{i} - min(depthAll{i}), 150, 'marker', '+', ...
            'linewidth', 3, 'markeredgecolor', [.3 .7 .7]);
        set(gca, 'fontsize', 14);
        xlabel('stress (psi)');
        ylabel('aspiration depth (\mum)');
        title(['E' num2str(embryoNums{i})]);
        grid on;
        xlim([0 .05]);
        ylim([0 20]);
        plotHandles = [plotHandles h];
        
%         p = polyfit(pressureAll{i}, depthAll{i}, 1);
%         depthFit = polyval(p, pressureAll{i});
        
%         hold on;
%         plot(pressureAll{i}, depthFit, 'linewidth', 2, 'color', 'k');
        
%         [fileName, pathName] = uigetfile([fileDir '\E' num2str(embryoNums{i}) '*.jpg']);
%         A = rgb2gray(imread([pathName fileName]));
%         figure(f6);
%         clf;
%         imshow(A);
%         movegui(f6, [-100, -100]);
%     
        % find pressure for Lp = Rp, and subtract "holding" pressure
        % eqPressure is in psi, cellSize is cell radius in microns when
        % cell is aspirated at equipibrium (Lp = Rp)
        dSub = depthAll{i} - min(depthAll{i});
        pSub = pressureAll{i} - min(pressureAll{i});
        
        pInv = polyfit(dSub, pSub, 1);
        pEq = polyval(pInv, 20); 
        eqPressure(i) = pEq;
        
        % units of Pa * m = N / m
        surfTension(i) = pEq*6895/(2*(1/(20*10^-6) - 1/(cellSize(i)*10^-6)));

    end
end

% xlim([-.05, .42]);
% ylim([0 .27]);
grid on;
% legend(plotHandles, {'E1', 'E2', 'E3', 'E4', 'E5'});

if isequal(type, 'Membrane')
    save([procDir '\' type '\all' type '.mat'], 'depthAll', 'pressureAll', ...
    'cellSize', 'surfTension', 'eqPressure');
else
    save([procDir '\' type '\all' type '.mat'], 'depthAll', 'pressureAll', ...
    'cellSize');
end




