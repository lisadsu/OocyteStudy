% analyze microinjection data from 9-24-14
% Look at embryo mechanics, not CG images
% Livia Zarnescu
% 1-27-15

%% 1. Load in data

%clear all;
%close all;
plotGroup = [1 0 0 0];

groupNum = 4*ones(1,104);
groupNum(1:34) = 1;
groupNum(35:65) = 2;
groupNum(66:84) = 3;
embryoNum = 1:104;
% 
% [num, txt, ~] = xlsread('data/embryoPN-9-24-14.xlsx');
% pnList = num(:,2)';

colorMat = zeros(104,3);
k0list = [];
k1list = [];
n0list = [];
taulist = [];
n1list = [];

dataDir = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse embryo\';

for i = 1:104
    
    colorMat = fillColorMat(colorMat, groupNum(i), 1, i);
    embryoString = num2str(i);
    
    % load data from embryo i
    testPath = [dataDir, '9-24-14', ' analysis\AutoMeasure\aspiration_data_', ...
        '9_24_14', '_E', embryoString, '.mat'];
    
    if ~exist(testPath, 'file')
        testPath = [dataDir, '9-24-14', ' analysis\aspiration_data_', ...
            '9_24_14', '_E', embryoString, '.mat'];
    end
    
    % load params
    if exist(testPath, 'file')
        load(testPath);
        aspiration_depth = aspiration_depth * 40 * 10^-6 / 108; % convert from pixels to meters
        k0list = [k0list k0];
        k1list = [k1list k1];
        n0list = [n0list n0];
        taulist = [taulist tau];
        n1list = [n1list n1];
    else
        k0list = [k0list NaN];
        k1list = [k1list NaN];
        n0list = [n0list NaN];
        taulist = [taulist NaN];
        n1list = [n1list NaN];
    end

end

numsToPlot = ~isnan(k1list);
gN = groupNum(numsToPlot);
eN = embryoNum(numsToPlot);
k1N = k1list(numsToPlot);
n1N = n1list(numsToPlot);
tN = taulist(numsToPlot);
k0N = k0list(numsToPlot);
n0N = n0list(numsToPlot);
% pN = pnList(numsToPlot);

% 2. Plot data
figure;
for j = 1:4
    
    gP = gN;
    eP = eN;
    k1P = k1N;
    k0P = k0N;
    tP = tN;
    n1P = n1N;
    n0P = n0N;
    cP = colorMat(numsToPlot,:);
    
    for i = 1:4
        if (i ~= j)
            cP = cP(gP ~= i, :);
            k1P = k1P(gP ~= i);
            k0P = k0P(gP ~= i);
            n1P = n1P(gP ~= i);
            tP = tP(gP ~= i);
            n0P = n0P(gP ~= i);
            eP = eP(gP ~= i);
            gP = gP(gP ~= i);
        end
    end
    
    subplot(4,1,j);
    h = scatter3(k1P, n1P, k0P, 80, cP, 'filled');
    hold on;
    
    set(h, 'Marker', 'o');
    set(gca, 'FontSize', 14);
    title('3D scatter plot of parameters');
    xlabel('k1 parameter');
    ylabel('n1 parameter');
    zlabel('k0 parameter');
    view(0,48);
    set(gca, 'xscale', 'linear');
    set(gca, 'yscale', 'log');
    set(gca, 'zscale', 'log');
    xlim([.09 .19]);
    
    b = num2str(eP');
    c = cellstr(b);
    dx = -0.002; dy = 0.01; dz = .005; % displacement so the text does not overlay the data points
    text(k1P+dx, n1P+dy, k0P+dz, c);
    
end


%% 3. Load ref params

[paramsOut, mOut] = getMouseParams('mouse embryo', [zeros(1,25) 1 1 0]);

cOut = zeros(length(mOut), 3);
cOut(mOut == 7,:) = repmat([.3 .8 .3], length(mOut(mOut == 7)), 1);
cOut(mOut == 6,:) = repmat([.3 .3 .8], length(mOut(mOut == 6)), 1);

% make smoothed kernel density plot
[fN xiN] = ksdensity(paramsOut(mOut == 6, 1));
[fV xiV] = ksdensity(paramsOut(mOut == 7, 1));


figure;
set(gca, 'FontSize', 14);
title('Embryo viability after IVF');
xlabel('k1 parameter');
ylabel('density');

hold on;
v1 = plot(xiV, fV, 'Color', [.3 .8 .3], 'linewidth', 2);
nv1 = plot(xiN, fN, 'Color', [.3 .3 .8], 'linewidth', 2);
h = area(xiV, fV, 'EdgeColor', [.3 .8 .3], 'FaceColor', [.3 .8 .3]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .3);
h = area(xiN, fN, 'EdgeColor', [.3 .3 .8], 'FaceColor', [.3 .3 .8]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .3);
xlim([.05 .2]);
legend([v1, nv1], 'Viable', 'Nonviable');
grid on;
box on;

%% Make bar graph of mechanics for all groups

AbK = k1N(gN == 1);
ControlK = k1N(gN == 3);
NonControlK = k1N(gN == 4);

figure;
h1 = bar(1, mean(AbK), 'facecolor', [0 0 .6]);
hold on;
e1 = errorbar(1,mean(AbK),std(AbK),'color', 'k', 'linewidth', 2);
h2 = bar(2, mean(ControlK), 'facecolor', [0 .6 .6]);
e2 = errorbar(2,mean(ControlK),std(ControlK),'color', 'k', 'linewidth', 2);
h3 = bar(3, mean(NonControlK), 'facecolor', [0 .6 0]);
e3 = errorbar(3,mean(NonControlK),std(NonControlK),'color', 'k', 'linewidth', 2);
% h4 = bar(4, mean(RefK), 'facecolor', [0 .6 0]);
% e4 = errorbar(4,mean(RefK),std(RefK),'color', 'k', 'linewidth', 2);
set(gca, 'xtick', [1 2 3])
set(gca, 'fontsize', 14);
set(gca, 'xticklabel', {'Antibody', 'Control Injected', ...
    'Control Noninjected'});
ylabel('k_1');
title('Embryo Stiffness');
xlim([0.5 3.5]);
ylim([.1 .2]);
grid on;


%% kernel density plot of antibody vs control



% make smoothed kernel density plot
[f1 x1] = ksdensity(k1N(gN == 1));
[f2 x2] = ksdensity(k1N(gN == 3));


figure;
set(gca, 'FontSize', 14);
title('Mechanical Properties After Fertilization');
xlabel('k1 parameter');
ylabel('histogram density');

hold on;
v1 = plot(x1, f1, 'Color', [.8 .6 .3], 'linewidth', 2);
nv1 = plot(x2, f2, 'Color', [.3 .3 .8], 'linewidth', 2);
h = area(x1, f1, 'EdgeColor', [.8 .6 .3], 'FaceColor', [.8 .6 .3]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .3);
h = area(x2, f2, 'EdgeColor', [.3 .3 .8], 'FaceColor', [.3 .3 .8]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .3);
plot([0 .25], [0 0], '--k');
xlim([.08 .25]);
grid on;
box on;

v = plot([.1353 .1353 .1553 .1553], [0 100 100 0], ...
    'color', [.3 .8 .3], 'linewidth', 2);
h = area([.1353 .1353 .1553 .1553], [0 100 100 0], ...
    'edgecolor', [.3 .8 .3], 'facecolor', [.3 .8 .3]);
hc = get(h, 'Children');
set(hc, 'facealpha', .3);
ylim([0 40]);

legend([v1, nv1 v], 'Antibody Injected', 'Control Injected', 'Viable');

% test for normality...1st both are normal according to lillietest
%[h p] = lillietest(k1N(gN == 1), .05, 'norm', 1e-3)
%[h p] = lillietest(k1N(gN == 3), .05, 'norm', 1e-3)

% test for difference in means of stiffness values
[h p] = ttest2(k1N(gN == 1), k1N(gN == 3))

% test to see if they come from different distributions
[h p] = kstest2(k1N(gN == 1), k1N(gN == 3))




