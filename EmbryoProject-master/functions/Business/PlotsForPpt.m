% make plot of stats on egg freezing
%
%

year = [2013 2014 2015 2016 2017];
nEggCryo = [75 85 105 151 202];
nPct = 100*[75/593 87/565 105/600 96/495 .3];
barWidth = .2;
colorToPlot = [0 0 1];


figure(1); clf; 
hold on;
set(gca, 'fontsize', 14);
xlabel('year');
% ylabel('percent of all cycles (%)');
title('Egg cryopreservation at LPCH REI');

paramToPlot = nEggCryo;


for i = 1:length(year)
    
    if i == 4
        patch('xdata', year(i) + [-1*barWidth barWidth barWidth -1*barWidth], ...
            'ydata', [0 0 paramToPlot(i) paramToPlot(i)], 'edgecolor', colorToPlot/2, ...
            'facecolor', 1/2 + 1/2*colorToPlot, 'facealpha', 1, 'edgealpha', 1, ...
            'linewidth', 2, 'linestyle', '--');
    elseif i > 4
        patch('xdata', year(i) + [-1*barWidth barWidth barWidth -1*barWidth], ...
            'ydata', [0 0 paramToPlot(i) paramToPlot(i)], 'edgecolor', colorToPlot/2, ...
            'facecolor', [1 1 1], 'facealpha', 1, 'edgealpha', 1, ...
            'linewidth', 2, 'linestyle', '--');
    else
        patch('xdata', year(i) + [-1*barWidth barWidth barWidth -1*barWidth], ...
            'ydata', [0 0 paramToPlot(i) paramToPlot(i)], 'edgecolor', colorToPlot/2, ...
            'facecolor',  1/2 + 1/2*colorToPlot, 'facealpha', 1, 'edgealpha', 1, ...
            'linewidth', 2);
    end
end

xlim([2012.5 max(year) + .5]);
% ylim([0 50]);
set(gca, 'xtick', year);
grid on;

