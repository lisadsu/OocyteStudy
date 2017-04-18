%% 1. load in data

eList = 1:10;
dataDir = 'C:\Users\Livia\Desktop\IVF\Raw Data\Videos\Human\HumanOocyteMeas\MeasuredData';

numList = [];
k1GV = [];
n1GV = [];
k1M = [];
n1M = [];
colorMatGV = [];
colorMatM = [];
currNum = 0;
legendParam = [];
cTextGV = {};
cTextM = {};
colorList = rand(10,3);

figure(1);
hold on;


for i = 1:length(eList)
   
    % load in GV and M for each oocyte
    
    currNum = currNum + 1;
    currColor = colorList(i,:);
    currFileName = ['E' num2str(eList(i)) '_GV'];
    
    if exist([dataDir '\' currFileName '.mat'], 'file')
        load([dataDir '\' currFileName '.mat']);
        k1GV = [k1GV k1];
        n1GV = [n1GV n1];
        cTextGV = {cTextGV{:}, ['E' num2str(i) '\_GV']};
        colorMatGV = [colorMatGV; currColor];
    end
    
    currFileName = ['E' num2str(eList(i)) '_M'];
    
    if exist([dataDir '\' currFileName '.mat'], 'file')
        load([dataDir '\' currFileName '.mat']);
        k1M = [k1M k1];
        n1M = [n1M n1];
        cTextM = {cTextM{:}, ['E' num2str(i) '\_M']};
        colorMatM = [colorMatM; currColor];
    end
            
end






%% 3. scatter plot


figure(3);
clf;
hold on;

scatter(k1GV, n1GV, 150, colorMatGV);
scatter(k1M, n1M, 150, colorMatM, 'filled');

for i = 1:length(eList)
   plot([k1GV(i), k1M(i)], [n1GV(i), n1M(i)], 'color', colorMatGV(i,:), 'linewidth', 2); 
end

set(gca, 'FontSize', 14);
title('human oocyte mechanics');
grid on;
xlabel('k_1 parameter');
ylabel('\eta_1 parameter');

dx = -0.01; dy = 0.015; % displacement so the text does not overlay the data points
% text(k1GV+dx, n1GV+dy, cTextGV'); 
% text(k1M+dx, n1M+dy, cTextM'); 

figure(4);
clf; hold on;
scatter(k1M - k1GV, n1M - n1GV, 150, colorMatGV, 'filled');
set(gca, 'FontSize', 14);
title('human oocyte mechanics');
grid on;
xlabel('k_1 parameter');
ylabel('\eta_1 parameter');