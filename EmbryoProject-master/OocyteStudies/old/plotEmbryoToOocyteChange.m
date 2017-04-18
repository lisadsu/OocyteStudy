% Plot transition from oocyte to embryo (before and after IVF)
% Livia Zarnescu
% 8/1/14

% clear all;
% close all;

%
whichToPlotEmbryo = [zeros(1,8) ones(1,5) 0 1 1 1 1 0 1 1 0 0 0 0 0 0 0 1];
whichToPlotOocyte = [0 0 1 1 1 1 1 0 1];

[embryoP, embryoM] = getMouseParams('mouse embryo', whichToPlotEmbryo);
[oocyteP, oocyteM] = getMouseParams('mouse oocyte', whichToPlotOocyte);

% set colors to plot
colorMatEmbryo = .1*ones(length(embryoM), 3);
colorMatEmbryo(embryoM == 2, 3) = .6; % bad are [.1 .1 .5]
colorMatEmbryo(embryoM == 4, 2) = .6; % good are [.1 .5 .1]

colorMatEmbryo(embryoM == 5, 1) = .8; % unknown are [.8 .8 .1]
colorMatEmbryo(embryoM == 5, 2) = .8;
colorMatEmbryo(embryoM == 6, 1) = .9; % bad are [.9, .1, .1]
colorMatEmbryo(embryoM == 7, 2) = .9; % good are [.1, .9, .1]

colorMatOocyte = .3*ones(length(oocyteM), 3);
colorMatOocyte(oocyteM == 1, 1) = .8; % GV are [.8 .5 .5]
colorMatOocyte(oocyteM == 2, 3) = .8; % bad M2 are [.5 .5 .8]
colorMatOocyte(oocyteM == 4, 2) = .8; % good M2 are [.5 .8 .5]
colorMatOocyte(oocyteM == 3, 2) = .8; % M1 are [.5 .8 .8]
colorMatOocyte(oocyteM == 3, 3) = .8;

colorMatOocyte(oocyteM == 5, 1) = .9; % unknown are [.9, .9, .5]
colorMatOocyte(oocyteM == 5, 2) = .9;
colorMatOocyte(oocyteM == 6, 1) = .9; % bad are [.9, .5, .5]
colorMatOocyte(oocyteM == 7, 2) = .9; % good are [.5, .9, .5]

figure;
hE = scatter3(embryoP(:,1), embryoP(:,2), embryoP(:,4), 80, ...
    colorMatEmbryo);%, 'filled');
hold on;
hO = scatter3(oocyteP(:,1), oocyteP(:,2), oocyteP(:,4), 80, ...
    colorMatOocyte);%, 'filled');

% plot trajectories
embryoUnkP = embryoP(embryoM >= 5, :);
embryoUnkM = embryoM(embryoM >= 5);
oocyteUnkP = oocyteP(oocyteM >= 5, :);
oocyteUnkM = oocyteM(oocyteM >= 5);
colorMatEmbryoUnk = colorMatEmbryo(embryoM >= 5, :);
colorMatOocyteUnk = colorMatOocyte(oocyteM >= 5, :);

if (length(embryoUnkM) ~= length(oocyteUnkM))
    ME = MException('plotEmbryoToOocyteChange:unequalNumbersToCompare', ...
        ['Must have equal numbers of oocytes and embryos with morphology' ...
        'value 5 to properly draw trajectories']);
    throw(ME);
end

% for i = 1:length(embryoUnkM)
%
%    hold on;
%    vectarrow(oocyteUnkP(i, [1 2 4]), embryoUnkP(i, [1 2 4]));
%
% end

line([oocyteUnkP(:,1), embryoUnkP(:,1)]', [oocyteUnkP(:,2), embryoUnkP(:,2)]', ...
    [oocyteUnkP(:,4), embryoUnkP(:,4)]', 'color', [0 0 0]);

scatter3(embryoUnkP(:,1), embryoUnkP(:,2), embryoUnkP(:,4), 80, colorMatEmbryoUnk, 'filled');
hold on;
scatter3(oocyteUnkP(:,1), oocyteUnkP(:,2), oocyteUnkP(:,4), 80, colorMatOocyteUnk, 'filled');
scatter3(embryoUnkP(:,1), embryoUnkP(:,2), embryoUnkP(:,4), 80, [0 0 0]);
set(gca, 'xscale', 'linear');
set(gca, 'yscale', 'log');
set(gca, 'zscale', 'log');
set(gca, 'fontsize', 14);

a = (1:53)';
a([5 31 34]) = NaN;
a = a(~isnan(a));
b = num2str(a);
c = cellstr(b);
dx = -0.004; dy = 1.1; dz = 1.1; % displacement so the text does not overlay the data points
he = text(embryoUnkP(:,1)+dx, embryoUnkP(:,2)*dy, embryoUnkP(:,4)*dz, c, 'fontsize', 14');
ho = text(oocyteUnkP(:,1)+dx, oocyteUnkP(:,2)*dy, oocyteUnkP(:,4)*dz, c, 'fontsize', 14');

ylim([.2 20]);
zlim([.005 .3]);

%% plot mean embryo-oocyte stiffness change

figure;
h1 = bar(1, mean(oocyteUnkP(:,1)), 'facecolor', [.8 .1 .1]);
hold on;
e1 = errorbar(1,mean(oocyteUnkP(:,1)),std(oocyteUnkP(:,1)),'color', 'k', 'linewidth', 2);
h2 = bar(2, mean(embryoUnkP(:,1)), 'facecolor', [.1 .1 .8]);
e2 = errorbar(2,mean(embryoUnkP(:,1)),std(embryoUnkP(:,1)),'color', 'k', 'linewidth', 2);
set(gca, 'xtick', [1 2])
set(gca, 'fontsize', 14);
ylim([0 .2]);
set(gca, 'xticklabel', {'Oocyte', 'Embryo'});
ylabel('Cell Stiffness');
title('Zona hardening during fertilization');
xlim([0.5 2.5]);
grid on;


%% Plot zona hardening vs starting stiffness

colorMat = zeros(length(oocyteUnkM), 3);
colorMat(oocyteUnkM == 5, :) = repmat([.8 .8 .3], sum(oocyteUnkM == 5), 1);
colorMat(oocyteUnkM == 6, :) = repmat([.8 .3 .3], sum(oocyteUnkM == 6), 1);
colorMat(oocyteUnkM == 7, :) = repmat([.3 .8 .3], sum(oocyteUnkM == 7), 1);

figure, scatter(oocyteUnkP(:,1), embryoUnkP(:,1)./oocyteUnkP(:,1), 100, ...
    [.7 .7 0], 'filled');
% ylim([0 3]);
set(gca, 'fontsize', 14);
xlabel('Stiffness before fertilization');
ylabel('Change in stiffness');
title('Degree of zona hardening during fertilization');
grid on;
xlim([.04 .11]);











