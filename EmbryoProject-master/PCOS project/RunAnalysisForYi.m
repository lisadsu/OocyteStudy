% Run analysis of follicle distribution on Yi's data from Aaron's lab
% LZ 8/5/15

% First, format files appropriately for input into RipleyGUI

dataFolder = 'C:\Users\Livia\Desktop\IVF\Raw Data\Yi\RawData';
inputFolder = 'C:\Users\Livia\Desktop\IVF\Raw Data\Yi\InputFiles';

fileNameInput = {'d10 primary', 'd10 primordial', 'd10 secondary', ...
    'd21 primary', 'd21 antral', 'd21 primordial', 'd21 secondary', ...
    'P preovulatory', 'P antral', 'P primary', 'P primordial', ...
    'P secondary'};

for i = 5%length(fileNameInput)
    
    display(['current file name: ' fileNameInput{i}]);
    
    fileToLoad = uigetfile([dataFolder '\*.*']);
    [numX, txtX, rawX] = xlsread([dataFolder '\' fileToLoad], -1);
    [numY, txtY, rawY] = xlsread([dataFolder '\' fileToLoad], -1);
    [numZ, txtZ, rawZ] = xlsread([dataFolder '\' fileToLoad], -1);
    
    numX = (numX - mean(numX))/264.58;
    numY = (numY - mean(numY))/264.58;
    numZ = (numZ - mean(numZ))/264.58;
    
    if (length(numX) == length(numY) && length(numX) == length(numZ))
        
        % write X,Y,Z coordinates to text file
        [fid, message] = fopen([inputFolder '\' fileNameInput{i} '.txt'], 'w');
        
        % go line by line
        for j = 1:length(numX)
            fprintf(fid,'%3.3f,%3.3f,%3.3f\r\n', numX(j), numY(j), numZ(j));
        end
        
        % close text file
        fclose(fid);
    end
end

%%
% 
h = gcf;
axesObjs = get(h, 'children');
dataObjs = get(axesObjs, 'children');

% secondary
x3 = get(dataObjs{5}(1), 'xdata')
y3 = get(dataObjs{5}(1), 'ydata')

% primary
x2 = get(dataObjs{5}(4), 'xdata')
y2 = get(dataObjs{5}(4), 'ydata')

% primordial
x1 = get(dataObjs{5}(7), 'xdata')
y1 = get(dataObjs{5}(7), 'ydata')


% % for d21 stage 124 only
% h = gcf;
% axesObjs = get(h, 'children');
% dataObjs = get(axesObjs, 'children');
% 
% % antral
% x4 = get(dataObjs{5}(1), 'xdata')
% y4 = get(dataObjs{5}(1), 'ydata')


% % for P stage 345 only
% h = gcf;
% axesObjs = get(h, 'children');
% dataObjs = get(axesObjs, 'children');
% 
% % preovulatory
% x5 = get(dataObjs{5}(1), 'xdata')
% y5 = get(dataObjs{5}(1), 'ydata')
% 
% % antral
% x4 = get(dataObjs{5}(4), 'xdata')
% y4 = get(dataObjs{5}(4), 'ydata')

%% Plot

figure(5);
clf;
set(gca, 'fontsize', 14);
plot(x1, y1, 's-', 'color', [0 1 0], 'linewidth', 4, ...
    'markerfacecolor', [0 0 0], 'markeredgecolor', 'none');
hold on;
plot(x2, y2, 's-', 'color', [1 0 1], 'linewidth', 4, ...
    'markerfacecolor', [0 0 0], 'markeredgecolor', 'none');

plot(x3, y3, 's-', 'color', [0 0 1], 'linewidth', 4, ...
    'markerfacecolor', [0 0 0], 'markeredgecolor', 'none');

plot(x4, y4, 's-', 'color', [1 .5 0], 'linewidth', 4, ...
    'markerfacecolor', [0 0 0], 'markeredgecolor', 'none');

% plot(x5, y5, 's-', 'color', [1 0 0], 'linewidth', 4, ...
%     'markerfacecolor', [0 0 0], 'markeredgecolor', 'none');


% legend('primordial', 'primary', 'secondary', 'antral', 'preovulatory', ...
%     'location', 'northwest');
legend('primordial', 'primary', 'secondary', 'antral', ...
    'location', 'northwest');
% legend('primordial', 'primary', 'secondary', ...
%     'location', 'northwest');
grid on;
xlabel('d (\mum)');
ylabel('K(t) - E[K(t)]');
title('Day 21');





