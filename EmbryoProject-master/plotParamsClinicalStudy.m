% plot params for all patients in clinical study
% Livia Zarnescu Yanez 11-10-15

params = loadDataClinicalStudy();
numParticipants = params.numParticipants;
participantIDs = params.participantIDs;
numEmbryos = params.numEmbryos;
outcomeInfo = params.outcomeInfo;
zygoteM = params.zygoteM;
d3M = params.d3M;
blastM = params.blastM;
ICSI = params.ICSI;
PGD = params.PGD;
gender = params.gender;
patientAge = params.patientAge;
patientBMI = params.patientBMI;
totalNumEmbryos = sum(numEmbryos);
procDataPath = 'C:\Users\Livia\Desktop\IVF\Processed Data\Human\';

% init param vectors
mList = [];
k0list = [];
k1list = [];
n0list = [];
n1list = [];
taulist = [];
gList = [];
colorMat = [];
icsiList = [];
zList = [];
aList = [];
ptList = [];
currE = 0;
dayThreeM1 = [];
dayThreeM2 = [];
blastStage = [];
blastScore = [];
ageList = [];
bmiList = [];
aText = [];
cText = {};

ptsToPlot = [ones(1,26) 1 1 1 1 1*ones(1,3)]; 


% load in params by participant and embryo num
for i = 1:numParticipants
    if ptsToPlot(i)
        i
        
        for j = 1:numEmbryos(i)

            j
            currE = currE + 1;
            currDataPath = [procDataPath participantIDs{i} '\' ...
                participantIDs{i} '_E' num2str(j) '.mat'];
            ptList = [ptList i];
            
            % save embryo params and color
            if exist(currDataPath, 'file')%  && (ICSI{i}(j) == 1) 
                load(currDataPath);
                mList = [mList outcomeInfo{i}(j)>0];
                k0list = [k0list k0];
                k1list = [k1list k1];
                n0list = [n0list n0];
                n1list = [n1list n1];
                taulist = [taulist tau];
                gList = [gList gender{i}(j)];
                icsiList = [icsiList ICSI{i}(j)];
                zList = [zList zygoteM{i}(j)];
                aText = [aText currE];
                cText = {cText{:}, ['P', num2str(i), '\_E', num2str(j)]};
%                 cText = {cText{:}, num2str(j)};
                
                aPad = padarray(A,[0,65-length(A)],'replicate', 'post');
                aList = [aList; aPad(1:65)];
                ageList = [ageList patientAge(i)];
                bmiList = [bmiList patientBMI(i)];
                1
                % get day 3 morphology data
                dayThreeM = d3M{i}{j};
                if numel(dayThreeM) == 0
                    dayThreeM1 = [dayThreeM1 0];
                    dayThreeM2 = [dayThreeM2 0];
                    k1list(end) = NaN;
                else
                    dayThreeM1 = [dayThreeM1 str2num(dayThreeM(1))];
                    dayThreeM2 = [dayThreeM2 getNumFromNumeral(dayThreeM(2:end))];
                end
                2
                % get day 5 morphology data
                blastString = blastM{i}{j};
                if numel(blastString) == 0
                    blastStage = [blastStage 0];
                    blastScore = [blastScore 0];
                else
                    blastStage = [blastStage str2num(blastString(1))];
                    if numel(blastString) < 3
                        blastScore = [blastScore 0];
                    else
                        blastScore = [blastScore (double((blastString(2)-65)+(blastString(3)-65))/2)];
                    end
                end
                
                
                if ((blastStage(end) > 0) && (blastScore(end) < 1.6)) 
                    colorMat = [colorMat; [.2 .6 .9]];
                elseif (blastStage(end) > 0)
%                     k1list(end) = NaN;
%                     n1list(end) = NaN;
                    colorMat = [colorMat; [.85 .65 .2]];%[.9 .2 .8]];%
                else
                    colorMat = [colorMat; [.85 .65 .2]];
                end
                
%                 if outcomeInfo{i}(j) == 2 % unknown
%                     mList(end) = 2;
%                     colorMat(end,:) = [.8 .2 .8];
%                 end
%                 
                % implantation/pregnancy color coding
                if outcomeInfo{i}(j) == 10 % no hCG rise
                    mList(end) = 2;
                    colorMat(end,:) = [1 .3 .3];
                end
                
                if outcomeInfo{i}(j) == 11 % hCG rise, no pregnancy
                    mList(end) = 3;
                    colorMat(end,:) = [.6 .2 1];
                end
                
                if outcomeInfo{i}(j) == 12 % hCG rise, pregnancy
                    mList(end) = 4;
                    colorMat(end,:) = [.3 1 .3];
                end

                    % ICSI / IVF color coding
%                     if ICSI{i}(j) == 1
%                         colorMat(end,:) = [.2 .6 .9];
%                     elseif ICSI{i}(j) == 0
%                         colorMat(end,:) = [.85 .65 .2];
%                     end
                    
                    % PGD color coding
%                   if PGD{i}(j) == 1
%                       colorMat(end,:) = [.2 .6 .9]; % euploid
%                   elseif PGD{i}(j) == 0
%                       colorMat(end,:) = [.85 .65 .2]; % aneuploid
%                   else 
%                       k1list(end) = NaN;
%                   end

%                 % age / bmi color coding
%                 if patientBMI(i) > 25
%                     colorMat = [colorMat; [0 0 1]];
%                 else
%                     colorMat = [colorMat; [1 0 0]];
%                 end


                % remove D3 transfer data for P10_E3
                if (i == 10) && (j == 3)
                    k1list(end) = NaN;
                end
                
            else
%                 colorMat = [colorMat; [NaN NaN NaN]];
%                 mList = [mList NaN];
%                 k0list = [k0list NaN];
%                 k1list = [k1list NaN];
%                 n0list = [n0list NaN];
%                 n1list = [n1list NaN];
%                 taulist = [taulist NaN];
%                 icsiList = [icsiList NaN];
%                 aList = [aList; NaN*zeros(1,65)];
%                 gList = [gList NaN];
%                 zList = [zList NaN];
%                 aText = [aText NaN];
%                 cText = {cText{:}, ''};
%                 ageList = [ageList NaN];
%                 bmiList = [bmiList NaN];
%                 dayThreeM1 = [dayThreeM1 NaN];
%                 dayThreeM2 = [dayThreeM2 NaN];
%                 blastStage = [blastStage NaN];
%                 blastScore = [blastScore NaN];
            end
        end
    end
end


aList = 10^6*aList;
% colorMat(icsiList == 1,:) = repmat([NaN NaN NaN], length(icsiList(icsiList == 1)), 1);

transfer1 = k1list(mList==10);
transfer2 = n1list(mList==10);
transferC = colorMat(mList==10,:);

% now plot!
p1 = k1list;
p2 = n1list; %n1list;
p3 = dayThreeM1; %zList + .1*randn(1,length(zList));

% new plot

figure(1);
clf;
hold on;
h = cell(1,length(p1));
hLegend = cell(1,5);
theta = 0:0.1:(2*pi);

for i = 1:length(p1)
    if ~isnan(p1(i)) && mList(i) < 2
        
        xPatch = cos(theta)/60 + p1(i);
        yPatch = sin(theta)/30 + p2(i);
        h = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', colorMat(i,:)/2, ...
            'facecolor', colorMat(i,:), 'facealpha', .6, 'edgealpha', .9);

%             [xPatch, yPatch, zPatch] = sphere(10);
%             xPatch = xPatch/40 + p1(i);
%             yPatch = yPatch/20 + p2(i);
%             zPatch = zPatch/2 + p3(i);
%             currPoint = surf2patch(xPatch, yPatch, zPatch);

%             h = patch(currPoint, 'edgecolor', colorMat(i,:)/2, ...
%                 'facecolor', colorMat(i,:), 'facealpha', faceAlpha, 'edgealpha', .1);
%         
        
        hold on;        
        if isempty(hLegend{mList(i)+1})
            hLegend{mList(i)+1} = h;
        end
    end
end

for i = 1:length(p1)
    if ~isnan(p1(i)) && mList(i) > 1

        xPatch = cos(theta)/60 + p1(i);
        yPatch = sin(theta)/30 + p2(i);
        h = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', colorMat(i,:)/2, ...
            'facecolor', colorMat(i,:), 'facealpha', 1, 'edgealpha', .9);
        hold on;        
        if isempty(hLegend{mList(i)+1})
            hLegend{mList(i)+1} = h;
        end
    end
end


set(gca, 'FontSize', 14);
grid on;
axis([.3 .8 0.05 .85]);
xlabel('k_1 parameter');
ylabel('\eta_1 parameter');

% legend([hLegend{1}, hLegend{2}], ...
%     {'Arrested', 'Blastocyst'}, 'Location', 'NorthEast');


dx = -0.01; dy = 0.015; dz = .5; % displacement so the text does not overlay the data points
hold on;
currPt = 29;
% text(p1+dx, p2+dy, cText'); %, p3+dz, c);
% text(p1(ptList == currPt)+dx, p2(ptList == currPt)+dy, {cText{ptList == currPt}}'); %, p3+dz, c);
% view(.5,90)


%% 1. Plot ICSI vs IVF mech properties

figure(2);
clf;
k11 = bar(1, mean(k1list(icsiList == 0 & ~isnan(k1list))), .8, 'facecolor', [0 .6 .6]);
hold on;
% ek11 = errorbar(1, mean(k1list(icsiList == 0)),std(k1list(icsiList == 0)),'color', 'k', 'linewidth', 2);
ek11 = terrorbar(1,mean(k1list(icsiList == 0 & ~isnan(k1list))), ...
    std(k1list(icsiList == 0 & ~isnan(k1list))), .1);
set(ek11, 'color', 'k', 'linewidth', 2);

k12 = bar(1.8, mean(k1list(icsiList == 1 & ~isnan(k1list))), .8, 'facecolor',[.6 .9 .9] );
ek12 = terrorbar(1.8, mean(k1list(icsiList == 1 & ~isnan(k1list))), ...
    std(k1list(icsiList == 1 & ~isnan(k1list))), .1);
set(ek12, 'color', 'k', 'linewidth', 2);

n11 = bar(3.2, mean(n1list(icsiList == 0 & ~isnan(k1list))), .8, 'facecolor', [1 .7 .2]);
en11 = terrorbar(3.2, mean(n1list(icsiList == 0 & ~isnan(k1list))), ...
    std(n1list(icsiList == 0 & ~isnan(k1list))), .1);
set(en11, 'color', 'k', 'linewidth', 2);

n12 = bar(4, mean(n1list(icsiList == 1 & ~isnan(k1list))), .8, 'facecolor',[1 .9 .6]);
en12 = terrorbar(4, mean(n1list(icsiList == 1 & ~isnan(k1list))), ...
    std(n1list(icsiList == 1 & ~isnan(k1list))), .1);
set(en12, 'color', 'k', 'linewidth', 2);

set(gca, 'fontsize', 14);
ylim([0 1]);
set(gca, 'xtick', [1.4 3.6]);
set(gca, 'xticklabel', {'k_1', '\eta_1'});
% xlabel('fertilization method');
ylabel('parameter value');
% title('fertilization method affects mechanics');
xlim([0.5 4.5]);

[p h] = ranksum(k1list(icsiList == 0), k1list(icsiList == 1))
[p h] = ranksum(n1list(icsiList == 0), n1list(icsiList == 1))

%% 2. Plot all aspiration depth curves

tVec = (0:1:64)/60;
figure(3);
clf;


for i = 1:length(mList)
    if ~isnan(mList(i))
        plot(tVec, 10^6*aList(i,:), 'color', colorMat(i,:), 'marker', 'o');
        hold on;
    end
end

axis([0 max(tVec) 0 45]);

%% 4. Extract and plot morphology



figure(1); clf;
h = scatter(dayThreeM1 + .75*rand(1,length(dayThreeM1)), ...
    dayThreeM2  + .75*rand(1,length(dayThreeM1)), 200, colorMat, 'filled');

set(h, 'Marker', 'o');
set(gca, 'FontSize', 14);
title('Mechanical parameters predict blastocyst formation');
xlabel('day3 morphology parameter 1');
ylabel('day3 morphology parameter 2');
% zlabel('zygote morphology');
% axis([min(p1) max(p1) min(p2) max(p2) min(p3) max(p3)]);
% set(gca, 'yscale', 'log')
% axis([min(p1) max(p1) min(p3) max(p3)]);
% set(gca, 'zscale', 'linear');
grid on;
% xlim([.3 .7]);
% ylim([0 .8]);

% b = num2str(aText');
% c = cellstr(b);
% dx = -0.01; dy = 0.015; dz = .5; % displacement so the text does not overlay the data points
% hold on;
% text(p1+dx, p2+dy,p3+dz, cText'); %, p3+dz, c);
% view(.5,90)


%% 5. Bar chart comparison of morphology vs mechanics

figure(4);
clf;
% day 3 morphology ROC
h1 = bar(.5, 0.656, .6, 'facecolor', [.8 .3 .3]);
hold on;
e1 = terrorbar(.5,0.656,.025,.1);
set(e1, 'color', 'k', 'linewidth', 2);

% mechanics ROC
h2 = bar(1.5, .811, .6, 'facecolor', [.3 .7 .8]);
e2 = terrorbar(1.5,.811,.011,.1);
set(e2, 'color', 'k', 'linewidth', 2);

% combined ROC
h3 = bar(2.5, .865, .6, 'facecolor', [.6 .3 .6]);
e3 = terrorbar(2.5,.865,.014,.1);
set(e3, 'color', 'k', 'linewidth', 2);

% % day 3 morphology PR
% h3 = bar(1.85, 0.599, .3, 'facecolor', [.8 .3 .3]);
% e3 = errorbar(1.85,0.599,.045);
% set(e1, 'color', 'k', 'linewidth', 2);
% 
% % mechanics PR
% h4 = bar(2.15, .703, .3, 'facecolor', [.3 .7 .8]);
% e4 = errorbar(2.15, .703, .016);
% set(e1, 'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [])
ylim([.5 1]);
set(gca, 'fontsize', 14);
% set(gca, 'xticklabel', {'ROC', 'PR'});
% ylabel('Area under curve');
% title('Comparison of mechanics and day 3 morphology');
xlim([0 3]);
grid on;

%% 6. SVM figure plot

figure(1);
% clf;
% h = scatter(p1,p2,200,colorMat,'linewidth',4);
hold on;

axisLims = [.3 .8 .05 .85];
embryoClassifier = fitcsvm([k1list; n1list]', mList == 1, 'ResponseName', ...
    'Viability', 'KernelFunction', 'rbf', 'KernelScale', 1, 'standardize', true);
[x1Grid,x2Grid] = meshgrid(linspace(axisLims(1),axisLims(2),100),...
    linspace(axisLims(3),axisLims(4),100));
xGrid = [x1Grid(:),x2Grid(:)];
[~,scores] = predict(embryoClassifier,xGrid);
[c, hh] = contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),-.1*[1 1],'k', 'linewidth', 2);

axis([.38 .81 0 .75]);

% set(h, 'Marker', 'o');
set(gca, 'FontSize', 14);
title('IVF');
xlabel('k1 parameter');
ylabel('n1 parameter');
grid on;
axis(axisLims);



%% Plot aspiration depth only

% load in aspiration depth data
load(['C:\Users\Livia\Desktop\IVF\Processed Data\Human' ...
    '\MECH029\MECH029_E6.mat'])
xdata = t;
ydata = A*10^6;

% plot
x1 = xdata(2:end);
x2 = [xfine - min(xfine)];
y1 = ydata(2:end);
y2 = [10^6*yfit];
x1 = [x1 NaN*ones(1, length(x2) - length(x1))];
y1 = [y1 NaN*ones(1, length(y2) - length(y1))];

figure(1); clf;
h = plot(x1, y1, 'color', [0 0 1], 'marker', 'o', 'linestyle', 'none');
hold on;
h = plot(x2, y2, 'color', [0 0 1], 'linestyle', '-');

axis([0 .65 0 24]);
set(gca, 'fontsize', 14);
grid on

%% Make heatmap of blast vs no blast

figure(4); clf;
spacing = .02;
xVals = .25:spacing:.8;
yVals = .05:spacing:.8;

valInRange = ((k1list > .3) & (k1list < .8) & (n1list > .05) & (n1list < .8));
kAll = k1list(~isnan(mList) & ~isnan(k1list) & valInRange);
nAll = n1list(~isnan(mList) & ~isnan(k1list) & valInRange);
mAll = mList(~isnan(mList) & ~isnan(k1list) & valInRange);

nPoints = 10;
minDist = 0.1;
scalingFactor = 100;
blastGrid = zeros(length(xVals)-1, length(yVals)-1);

% make color mat
colorMapMat = zeros(100,3);
minColor = [.85 .65 .2];
maxColor = [.2 .2 .9];

for i = 1:100
   colorMapMat(i,:) = minColor*(100-i)/100 + maxColor*i/100;
end


% compute mean viability of up to nPoints nearest neighbors within minDist
for i = 1:size(blastGrid,1)
    for j = 1:size(blastGrid,2)
        
        % take centroid of each grid block
        kCurr = mean([xVals(i), xVals(i+1)]);
        nCurr = mean([yVals(j), yVals(j+1)]);
        
        % compute distance to all data points
        distAll = sqrt((kAll - kCurr).^2 + (nAll - nCurr).^2);        
        
        % kNN
%         [distSort, I] = sort(distAll, 2, 'ascend');
%         mSort = mAll(I);
%         
%         % only keep the first nPoints which are within minDist
%         distClose = distSort(1:nPoints);
%         mClose = mSort(1:nPoints);
%         mClose = mClose(distClose < minDist);
        
%         if length(mClose) < 2
%             meanViability = 0;
%         else
%             meanViability = mean(mClose);
%         end
% 
%         blastGrid(i,j) = meanViability;

        % gaussian kernel regression
        mScale = mAll*2 - 1;
        
        % weighted sum based on distance
        meanViability = sum(mAll .* exp(-1*(distAll.^2)*scalingFactor) ./...
            repmat(sum(exp(-1*(distAll.^2)*scalingFactor)), 1, length(mAll)));
        
        blastGrid(i,j) = meanViability;
        currColor = minColor*(1-meanViability) + maxColor*meanViability;
        
        xPatch = [xVals(i) xVals(i+1) xVals(i+1) xVals(i) xVals(i)];
        yPatch = [yVals(j) yVals(j) yVals(j+1) yVals(j+1) yVals(j)];
        h = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', currColor/2, ...
            'facecolor', currColor, 'facealpha', .8, 'edgealpha', .2);

         
        
    end
end

% blastGridSmooth = imfilter(blastGrid, fspecial('gaussian', [5 5], 2));
% blastGridSmooth = imadjust(blastGridSmooth);
% 
% 
% for i = 1:size(blastGridSmooth,1)
%     for j = 1:size(blastGridSmooth,2)    
%     
%         meanViability = blastGridSmooth(i,j);
%         currColor = minColor*(1-meanViability) + maxColor*meanViability;
% 
%         xPatch = [xVals(i) xVals(i+1) xVals(i+1) xVals(i) xVals(i)];
%         yPatch = [yVals(j) yVals(j) yVals(j+1) yVals(j+1) yVals(j)];
%         h = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', currColor/2, ...
%             'facecolor', currColor, 'facealpha', .8, 'edgealpha', .2);
%         
%         
%     end
% end

% xIncrement = xVals(2) - xVals(1);
% yIncrement = yVals(2) - yVals(1);
% [xGrid, yGrid] = meshgrid(xVals(1:end-1) + xIncrement, ...
%     yVals(1:end-1) + yIncrement);


% colormap(colorMapMat);
% 
% % now plot grid
% surf(xGrid, yGrid, blastGrid');
% view(0,90);
axis([.25 .8 .05 .75]);

%% Plot ROC curves for mech vs morphology

figure(1); clf;
set(gca, 'FontSize', 14);
hold on;

% mech
h1 = plot([X1 1], [Y1 1], 'Color', [.3 .7 .8], 'LineWidth', 2);

% day 3 morphology
h2 = plot([X2 1], [Y2 1], 'Color', [.8 .3 .3], 'LineWidth', 2);

% optimal combination
h3 = plot([X3 1], [Y3 1], 'Color', [.6 .3 .6], 'LineWidth', 2);

plot([0 1], [0 1], 'color', .3*ones(1,3), 'linestyle', '--', 'linewidth', 2);

xlabel('1 - Specificity'); ylabel('Sensitivity')
title('');
grid on;
axis([0 1 0 1]);
legend([h1 h2 h3], {'Mechanics', 'Morphology', 'Both'})

save('clinicalROCdata.mat', 'X1', 'X2', 'X3', 'Y1', 'Y2', 'Y3', ...
    'AUC1', 'AUC2', 'AUC3');

% Mech only: sens 79%, spec 70%, ppv 72%, acc 74%, AUC 0.81
% Morph only: sens 68%, spec 59%, ppv 61%, acc 63%, AUC 0.66
% both: sens 77%, spec 79%, ppv 78%, acc 78%, AUC 0.86



%% Mech feature selection plot



figure(1); clf;
set(gca, 'FontSize', 14);
hold on;

plot([1 2 3 4], [AUC1, AUC2, AUC3, AUC4], 'color', [0 0 1], 'linewidth', 2);


xlabel('number of parameters'); ylabel('AUC_R_O_C')
title('');
grid on;
axis([1 4 .5 1]);






