% system identification applied to mouse embryo aspiration curves
% Applies 5 different models to the data:
% Maxwell (2 elements), Kelvin-Voigt (2 elements), SLS (3 elements),
% modified SLS (our model, 4 elements), and Wiechert (5 elements)
%
% Fit each model to each embryo, and calculate its error compared to all
% other embryos
%
% First: run the script plotAllMouseEmbryosSoFar
% This will result in 2 matrices aspirationDataN and aspirationDataV
% containing the data for all embryos

fileToSave = 'functions\TestModels\fitDataViableOnly.mat';
aspirationDataAll = aspirationDataV; %[aspirationDataV; aspirationDataN];
s = size(aspirationDataAll,1);
F = .325 * 6894.75729 * pi * (20 * 10^-6)^2; % pressure * area
plotInput = 1;

% training / in-sample error
fMaxwell = zeros(1,s);
fKelvinVoigt = zeros(1,s);
fSLS = zeros(1,s);
fMod = zeros(1,s);
fFive = zeros(1,s);

cvFolds = 10;

% testing / validation error (for each embryo, with cross validation)
vMaxwell = zeros(1,s);
vKelvinVoigt = zeros(1,s);
vSLS = zeros(1,s);
vMod = zeros(1,s);
vFive = zeros(1,s);


%% 1. Maxwell model

plotInput = 1;

for i = 1:s
    
    i 
    xData = aspirationDataAll(i,:);
    tData = 0:.013:(.013*(length(xData) - 1));

    startParams = zeros(1,2);
    startParams(1) = .15; % k
    startParams(2) = .5; % n
    
    ind = crossvalind('kfold', length(xData), cvFolds);
    xOutTraining = zeros(1,length(xData));
    xOutValidation = zeros(1,length(xData));
    
    for j = 0:cvFolds
    
        % in-sample
        if j == 0
            xTrain = xData;
            tTrain = tData;
            xTest = xData;
            tTest = tData;
        % CV
        else
            xTrain = xData(ind ~= j);
            tTrain = tData(ind ~= j);
            xTest = xData(ind == j);
            tTest = tData(ind == j);
        end
        
        % get best params for fit
        [paramsOut, fval] = MaxwellOptimizeParams(tTrain, xTrain, startParams, F, 0);
        k = paramsOut(1);
        n = paramsOut(2);
        params = [F, n];
        options = odeset();

        % evaluate ode at testing data points
        [~, xOut] = ode45( @(t,x)MaxwellTest(t,x,params), tTest, ...
            F/k, options);

        % fval is the in-sample fvalue
        if j == 0
            xOutTraining = xOut';
        else
        % accumulate test predictions
            xOutValidation(ind == j) = xOut';
        end
        
    end
    
    % place extra weight on good fit in first data points
    weighting = 1*ones(1, length(xData));

    if length(xData) > 20
        weighting(1:10) = 2;
        weighting(1) = 10;
    else
        weighting(1:length(xData)-2) = 1;
    end

    % calculate out-of-sample fvalue
%     fMaxwell(i) = sum(weighting.*(xOutTraining - xData).^2);
%     vMaxwell(i) = sum(weighting.*(xOutValidation - xData).^2);
    
    if (plotInput && i == 3)
        figure(1);
        clf;
        plot(tData, 10^6*xData, 'Marker', 'o', 'LineStyle', 'none');
        hold on;
        plot(tData, 10^6*xOutTraining);
        xlabel('time (seconds)');
        ylabel('aspiration depth (\mum)');
        ylim([10 20]);
        xlim([0 .5]);
    end
    
    
end

%% 2. Kelvin-Voigt model

plotInput = 1;

for i = 1:s
    
    i 
    xData = aspirationDataAll(i,:);
    tData = 0:.013:(.013*(length(xData) - 1));

    startParams = zeros(1,2);
    startParams(1) = .15; % k
    startParams(2) = .002; % n
    
    ind = crossvalind('kfold', length(xData), cvFolds);
    xOutTraining = zeros(1,length(xData));
    xOutValidation = zeros(1,length(xData));
    
    for j = 0:cvFolds
    
        % in-sample
        if j == 0
            xTrain = xData;
            tTrain = tData;
            xTest = xData;
            tTest = tData;
        % CV
        else
            xTrain = xData(ind ~= j);
            tTrain = tData(ind ~= j);
            xTest = xData(ind == j);
            tTest = tData(ind == j);
        end
        
        % get best params for fit
        [paramsOut, fval] = KelvinVoigtOptimizeParams(tTrain, xTrain, startParams, F, 0);
        k = paramsOut(1);
        n = paramsOut(2);
        params = [F, k, n];
        options = odeset();
        
        % evaluate ode at testing data points
        [~, xOut] = ode45( @(t,x)KelvinVoigtTest(t,x,params), tTest, ...
            0, options);

        % fval is the in-sample fvalue
        if j == 0
            xOutTraining = xOut';
        else
        % accumulate test predictions
            xOutValidation(ind == j) = xOut';
        end
        
    end
    
    % place extra weight on good fit in first data points
    weighting = 1*ones(1, length(xData));

    if length(xData) > 20
        weighting(1:10) = 2;
        weighting(1) = 10;
    else
        weighting(1:length(xData)-2) = 1;
    end

    % calculate out-of-sample fvalue
%     fKelvinVoigt(i) = sum(weighting.*(xOutTraining - xData).^2);
%     vKelvinVoigt(i) = sum(weighting.*(xOutValidation - xData).^2);
    
    if (plotInput && i == 3)
        figure(2);
        clf;
        plot(tData, 10^6*xData, 'Marker', 'o', 'LineStyle', 'none');
        hold on;
        plot(tData, 10^6*xOutTraining);
        xlabel('time (seconds)');
        ylabel('aspiration depth (\mum)');
        ylim([10 20]);
        xlim([0 .5]);
    end
    
end


%% 3. SLS model


plotInput = 0;

for i = 1:s
    
    i 
    xData = aspirationDataAll(i,:);
    tData = 0:.013:(.013*(length(xData) - 1));

    startParams = zeros(3,1);
    startParams(1) = .05; % k0
    startParams(2) = .15; % k1
    startParams(3) = .005; % n
    
    ind = crossvalind('kfold', length(xData), cvFolds);
    xOutTraining = zeros(1,length(xData));
    xOutValidation = zeros(1,length(xData));
    
    for j = 0:cvFolds
    
        % in-sample
        if j == 0
            xTrain = xData;
            tTrain = tData;
            xTest = xData;
            tTest = tData;
        % CV
        else
            xTrain = xData(ind ~= j);
            tTrain = tData(ind ~= j);
            xTest = xData(ind == j);
            tTest = tData(ind == j);
        end
        
        % get best params for fit
        [paramsOut, fval] = SLSOptimizeParams(tTrain, xTrain, startParams, F, 0);
        k0 = paramsOut(1);
        k1 = paramsOut(2);
        n = paramsOut(3);
        params = [F, k0, k1, n];
        
        
        % evaluate ode at testing data points
        [~, xOut] = ode45( @(t,x)SLSTest(t,x,params), tTest, ...
            F/(k0 + k1), options);

        % fval is the in-sample fvalue
        if j == 0
            xOutTraining = xOut';
        else
        % accumulate test predictions
            xOutValidation(ind == j) = xOut';
        end
        
    end
    
    % place extra weight on good fit in first data points
    weighting = 1*ones(1, length(xData));

    if length(xData) > 20
        weighting(1:10) = 2;
        weighting(1) = 10;
    else
        weighting(1:length(xData)-2) = 1;
    end

    % calculate out-of-sample fvalue
    fSLS(i) = sum(weighting.*(xOutTraining - xData).^2);
    vSLS(i) = sum(weighting.*(xOutValidation - xData).^2);
    
    if (plotInput && i == 3)
        figure(3);
        clf;
        plot(tData, 10^6*xData, 'Marker', 'o', 'LineStyle', 'none');
        hold on;
        plot(tData, 10^6*xOutTraining);
        xlabel('time (seconds)');
        ylabel('aspiration depth (\mum)');
        ylim([10 20]);
        xlim([0 .5]);
    end
    
end



% 4. Modified SLS (our model)


plotInput = 0;

for i = 1:s
    
    i 
    xData = aspirationDataAll(i,:);
    tData = 0:.013:(.013*(length(xData) - 1));

    startParams = zeros(4,1);
    startParams(1) = .02; % k0
    startParams(2) = .15; % k1'
    startParams(3) = .001; % n0
    startParams(4) = 10; % n1
    
    ind = crossvalind('kfold', length(xData), cvFolds);
    xOutTraining = zeros(1,length(xData));
    xOutValidation = zeros(1,length(xData));
    
    for j = 0:cvFolds
    
        % in-sample
        if j == 0
            xTrain = xData;
            tTrain = tData;
            xTest = xData;
            tTest = tData;
        % CV
        else
            xTrain = xData(ind ~= j);
            tTrain = tData(ind ~= j);
            xTest = xData(ind == j);
            tTest = tData(ind == j);
        end
        
        % get best params for fit
        [paramsOut, fval] = ModifiedSLSOptimizeParams(tTrain, xTrain, startParams, F, 0);
        k0 = paramsOut(1);
        k1 = paramsOut(2);
        n0 = paramsOut(3);
        n1 = paramsOut(4);

        params = [F, k0, k1, n0, n1];
        options = odeset();
        x0 = F/(k0 + k1);
        
        % evaluate ode at testing data points
        [~, xOut] = ode45( @(t,x)ModifiedSLSTest(t,x,params), tTest, ...
            [0 x0 0 x0], options);

        % fval is the in-sample fvalue
        if j == 0
            xOutTraining = xOut(:,4)';
        else
        % accumulate test predictions
            xOutValidation(ind == j) = xOut(:,4)';
        end
        
    end
    
    % place extra weight on good fit in first data points
    weighting = 1*ones(1, length(xData));

    if length(xData) > 20
        weighting(1:10) = 2;
        weighting(1) = 10;
    else
        weighting(1:length(xData)-2) = 1;
    end

    % calculate out-of-sample fvalue
    fMod(i) = sum(weighting.*(xOutTraining - xData).^2);
    vMod(i) = sum(weighting.*(xOutValidation - xData).^2);
    
    if (i == 3 && plotInput)
        figure(4);
        clf;
        plot(tData, 10^6*xData, 'Marker', 'o', 'LineStyle', 'none');
        hold on;
        plot(tData, 10^6*xOutTraining);
        xlabel('time (seconds)');
        ylabel('aspiration depth (\mum)');
        ylim([10 20]);
        xlim([0 .5]);
    end
    
end


% 5. Five-param model (Wiechert)


plotInput = 0;

for i = 1:s
    
    i 
    xData = aspirationDataAll(i,:);
    tData = 0:.013:(.013*(length(xData) - 1));

    startParams = zeros(5,1);
    startParams(1) = .15; % k0
    startParams(2) = .03; % k1
    startParams(3) = .001; % k2
    startParams(4) = .5; % n0
    startParams(5) = .001; % n1
    
    ind = crossvalind('kfold', length(xData), cvFolds);
    xOutTraining = zeros(1,length(xData));
    xOutValidation = zeros(1,length(xData));
    
    for j = 0:cvFolds
    
        % in-sample
        if j == 0
            xTrain = xData;
            tTrain = tData;
            xTest = xData;
            tTest = tData;
        % CV
        else
            xTrain = xData(ind ~= j);
            tTrain = tData(ind ~= j);
            xTest = xData(ind == j);
            tTest = tData(ind == j);
        end
        
        % get best params for fit
        [paramsOut, fval] = FiveParamOptimizeParams(tTrain, xTrain, startParams, F, 0);
        k0 = paramsOut(1);
        k1 = paramsOut(2);
        k2 = paramsOut(3);
        n0 = paramsOut(4);
        n1 = paramsOut(5);

        params = [F, k0, k1, k2, n0, n1];
        options = odeset();
        x0 = F/(k0 + k1 + k2);

        % evaluate ode at testing data points
        [~, xOut] = ode45( @(t,x)FiveParamTest(t,x,params), tTest, ...
        [x0 0 x0 0 x0], options);

        % fval is the in-sample fvalue
        if j == 0
            xOutTraining = xOut(:,5)';
        else
        % accumulate test predictions
            xOutValidation(ind == j) = xOut(:,5)';
        end
        
    end
    
    % place extra weight on good fit in first data points
    weighting = 1*ones(1, length(xData));

    if length(xData) > 20
        weighting(1:10) = 2;
        weighting(1) = 10;
    else
        weighting(1:length(xData)-2) = 1;
    end

    % calculate out-of-sample fvalue
    fFive(i) = sum(weighting.*(xOutTraining - xData).^2);
    vFive(i) = sum(weighting.*(xOutValidation - xData).^2);
    
    if (i == 3 && plotInput)
        figure(5);
        clf;
        plot(tData, 10^6*xData, 'Marker', 'o', 'LineStyle', 'none');
        hold on;
        plot(tData, 10^6*xOutTraining);
        xlabel('time (seconds)');
        ylabel('aspiration depth (\mum)');
        ylim([10 20]);
        xlim([0 .5]);
    end
    
end


%% Evaluate fitting error

mean(fMaxwell)
mean(fKelvinVoigt)
mean(fSLS)
mean(fMod)
mean(fFive)


figure(6);
clf;
h1 = scatter([2 2 3 4 5], log([mean(fMaxwell), mean(fKelvinVoigt), ...
    mean(fSLS), mean(fMod), mean(fFive)]), 150, [0 0 1]);
set(h1, 'Marker', 'd');


% Evaluate validation error
vMeans = zeros(1,5);
vMeans(1) = mean(vMaxwell);
vMeans(2) = mean(vKelvinVoigt);
vMeans(3) = mean(vSLS);
vMeans(4) = mean(vMod);
vMeans(5) = mean(vFive);

figure(6);
hold on;
h2 = scatter([2 2 3 4 5], log(vMeans), 150, [1 0 0]);
set(h2, 'Marker', 'o');
set(gca, 'XTick', [2 3 4 5]);
set(gca, 'XTickLabel', {'2', '3', '4', '5'});
set(gca, 'fontsize', 14);
title('Comparison of Different Bulk Models');
xlabel('number of model parameters');
ylabel('log(error)');
xlim([1.9 5.1]);
ylim([-30 -16])
grid on;
legend([h1 h2], {'Training error', 'Testing error'})

% save(fileToSave, 'fMaxwell', 'fKelvinVoigt', 'fSLS', 'fMod', 'fFive', 'vMaxwell', 'vKelvinVoigt', 'vSLS', 'vMod', 'vFive');

[h p] = ranksum(vMod, vSLS)
[h p] = ranksum(vMod, vFive)












