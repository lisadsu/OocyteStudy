% temp calibration for ovary probe

% type K
actualTemp = [41.3,38.8,37.3,34.3,31.1,30.8,30.2,29.6,29.0,28.3,22.2];
measVoltage = [407,384,368,337,304,301,297,290,282,274,214];

figure(1);
clf;
set(gca, 'fontsize', 14);
hold on;
xlabel('voltage (mV)');
ylabel('temperature (C)');
scatter(measVoltage,actualTemp,100);

[p s] = polyfit(measVoltage,actualTemp,1);
tFit = polyval(p,measVoltage);
plot(measVoltage,tFit, 'k', 'linewidth', 1);