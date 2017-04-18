% analyze blast data

%% 1. load in data

blastList = 1:12;
repList = [3 4 3 5 4 4 5 3 5 4 4 4];
dataDir = 'C:\Users\Livia\Desktop\IVF\Raw Data\BlastMeasurements\MeasuredData';

numList = [];
k1list = [];
n1list = [];
colorMat = [];
currNum = 0;
legendParam = [];
cText = {};
colorList = rand(12,3);
colorList(11,:) = [1 0 0];

for i = 1:length(blastList)
    for j = 1:repList(i)
   
        currNum = currNum + 1;
        currFileName = ['B' num2str(blastList(i)) '_loc' num2str(j)];
        
        if exist([dataDir '\' currFileName '.mat'], 'file')
            load([dataDir '\' currFileName '.mat']);
        
            k1list = [k1list k1];
            n1list = [n1list n1];
            numList = [numList i];
            cText = {cText{:}, ['B' num2str(i) '\_L', num2str(j)]};
            
            legendParam = [legendParam i];
            colorMat = [colorMat; colorList(i,:)];
           
            
        end
            
    end
end




%% 2. scatter plot

p1 = k1list;
p2 = n1list;


figure(3);
clf;
hold on;
h = cell(1,length(k1list));
hLegend = cell(1,12);

colorMat(11,:) = [1 0 0];

for i = 1:length(k1list)
    i
    if ~isnan(legendParam(i))
        i
%         h{i} = plot(p1(i), p2(i), 'marker', 'o', 'markeredgecolor', ...
%             colorMat(i,:), 'markersize', 12, 'color', 'none', ...
%             'linewidth', 4);

        theta = 0:0.1:(2*pi);     
        xPatch = cos(theta)*(max(p1) - min(p1))/30 + p1(i);
        yPatch = sin(theta)*(max(p2) - min(p2))/25 + p2(i);
        

        h{i} = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', colorMat(i,:), ...
            'facecolor', colorMat(i,:), 'facealpha', .7);

        hold on;
        if isempty(hLegend{legendParam(i)})
            i
            hLegend{legendParam(i)} = h{i};
        end
    end
end

set(gca, 'FontSize', 14);
title('trophectoderm mechanics');
grid on;
% axis([.06 .145 -0.5 11]);
xlabel('k_1 parameter');
ylabel('\eta_1 parameter');
legend([hLegend{1}, hLegend{2}, hLegend{3}, hLegend{4}, hLegend{5}], ...
    {'B1', 'B2', 'B3', 'B4', 'B5'}, 'Location', 'NorthEast');


b = num2str(aText');
c = cellstr(b);
dx = -0.01; dy = 0.015; % displacement so the text does not overlay the data points
hold on;
% text(p1+dx, p2+dy, cText'); 

%% 4. consolidate by blast / morphology

kBlast = zeros(1,12);
nBlast = zeros(1,12);
mBlast = [1 1 0 0 1 0 0 0 1 1 1 1];
currColor = zeros(12,3);
mList = zeros(1, length(k1list));

% morphology order of 1-5 is 2,5,1,4,3
% morphology order of 8-12 is 9,12,11,10,8
% morphology order of 6/7 is ??
% I will give score 1 to top 3 from each 1-5 / 8-12 group

colorList = [.85 .65 .2; .2 .2 .6; .8 .2 .8];
cText = {};
figure(3);
clf;
theta = 0:.1:2*pi;
h = cell(1,2);

for i = 1:12
    kBlast(i) = median(k1list(legendParam == i));
    nBlast(i) = median(n1list(legendParam == i));
    currColor(i,:) = colorList(mBlast(i)+1,:);
    cText = {cText{:}, ['B' num2str(i)]};
    mList(legendParam == i) = repmat(mBlast(i), 1, sum(legendParam == i));
    
    xPatch = cos(theta)*.0022 + kBlast(i);
    yPatch = sin(theta)*1 + nBlast(i);
        
    h{mBlast(i)+1} = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', currColor(i,:), ...
            'facecolor',currColor(i,:), 'facealpha', .5, 'edgealpha', .5);
end


set(gca, 'fontsize', 14);
% scatter(kBlast, nBlast, 200, currColor, 'filled');
legend([h{1}, h{2}], {'poor', 'good'});
grid on;
xlabel('k_1 parameter');
ylabel('\eta_1 parameter');

dx = -0.002; dy = 0.5; % displacement so the text does not overlay the data points
hold on;
% text(kBlast+dx, nBlast+dy, cText'); 


%% 5. bar plot

% uses mBlast, kBlast and nBlast from previous section


figure(2);
clf;
paramToPlot = kBlast;
outcomeMeas = mBlast;

% "medium" 
k00 = bar(.5, mean(paramToPlot(outcomeMeas == 0)), .8, 'facecolor', colorList(1,:));
hold on;
ek00 = terrorbar(.5, mean(paramToPlot(outcomeMeas == 0)), ...
    std(paramToPlot(outcomeMeas == 0)), .1);
set(ek00, 'color', 'k', 'linewidth', 2);
% "good"
k01 = bar(1.5, mean(paramToPlot(outcomeMeas == 1)), .8, 'facecolor', colorList(2,:));
hold on;
ek01 = terrorbar(1.5, mean(paramToPlot(outcomeMeas == 1)), ...
    std(paramToPlot(outcomeMeas == 1)), .1);
set(ek01, 'color', 'k', 'linewidth', 2);



set(gca, 'xtick', [0.5 1.5])
set(gca, 'fontsize', 14);
ylim([0 .1]);
set(gca, 'xticklabel', {'poor', 'good'});
ylabel('k_1 parameter');
title('p = .03');
xlim([0 2]);
grid on;


[p h] = ranksum(paramToPlot(outcomeMeas == 0), paramToPlot(outcomeMeas == 1))


















