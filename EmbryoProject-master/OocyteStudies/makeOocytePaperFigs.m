% make figures for oocyte paper
%
% 


%% Figure 1. Sample aspiration curves for different stiffnesses

tList = [];
aList = [];
xList = [];
yList = [];

stiffColor = [0 .2 .8];
medColor = [.4 .4 .4];
softColor = [.8 .2 .8];

figure(1);
clf;
hold on;

% stiff
load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\1-20-16 analysis\aspiration_data_1_20_16_E14.mat')
tList = [tList; t(1:38)];
aList = [aList; A(1:38)];
xList = [xList; xfine];
yList = [yList; yfit];
plot(t(2:end), 10^6*A(1:end-1), 'marker', 'o', 'linestyle', 'none', 'color', stiffColor);
plot(xfine, 10^6*yfit, 'color', stiffColor);

load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\10-21-15 analysis\aspiration_data_10_21_15_E69.mat')
tList = [tList; t(1:38)];
aList = [aList; A(1:38)];
xList = [xList; xfine];
yList = [yList; yfit];
plot(t(2:end), 10^6*A(1:end-1), 'marker', 'o', 'linestyle', 'none', 'color', stiffColor);
plot(xfine, 10^6*yfit, 'color', stiffColor);

% med
load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\3-28-16 analysis\aspiration_data_3_28_16_E6.mat')
tList = [tList; t(1:38)];
aList = [aList; A(1:38)];
xList = [xList; xfine];
yList = [yList; yfit];
plot(t(2:end), 10^6*A(1:end-1), 'marker', 'o', 'linestyle', 'none', 'color', medColor);
plot(xfine, 10^6*yfit, 'color', medColor);

load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\1-20-16 analysis\aspiration_data_1_20_16_E9.mat')
tList = [tList; t(1:38)];
aList = [aList; A(1:38)];
xList = [xList; xfine];
yList = [yList; yfit];
plot(t(2:end), 10^6*A(1:end-1), 'marker', 'o', 'linestyle', 'none', 'color', medColor);
plot(xfine, 10^6*yfit, 'color', medColor);

% soft
load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\6-6-16 analysis\aspiration_data_6_6_16_E21.mat')
tList = [tList; t(1:38)];
aList = [aList; A(1:38)];
xList = [xList; xfine];
yList = [yList; yfit];
plot(t(2:end), 10^6*A(1:end-1), 'marker', 'o', 'linestyle', 'none', 'color', softColor);
plot(xfine, 10^6*yfit, 'color', softColor);

load('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\1-20-16 analysis\aspiration_data_1_20_16_E20.mat')
tList = [tList; t(1:38)];
aList = [aList; A(1:38)];
xList = [xList; xfine];
yList = [yList; yfit];
plot(t(2:end), 10^6*A(1:end-1), 'marker', 'o', 'linestyle', 'none', 'color', softColor);
plot(xfine, 10^6*yfit, 'color', softColor);


xlim([tList(1,2) 0.5]);
ylim([10 30])

grid on;
set(gca, 'fontsize', 14)
xlabel('time (s)');
ylabel('aspiration depth (\mum)');
title('Oocyte response to pressure');
