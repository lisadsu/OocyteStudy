% solve for system parameters using fminsearch
% 
% startParams = [n k]

function [paramsOut, fval] = KelvinVoigtOptimizeParams(tData, xData, startParams, F, plotInput)


[paramsOut, fval, exitflag] = fminsearch(@(paramsIn)KelvinVoigtTestFit(tData, ...
    xData, F, paramsIn), startParams, ...
    optimset('TolFun', 10^(-9), 'TolX', 10^(-5), 'MaxFunEvals', 1000));

if (plotInput)
    
    k = paramsOut(1);
    n = paramsOut(2);

    params = [F, k, n];
    options = odeset();

    % x(0) = 0;
    [tOut, xOut] = ode45( @(t,x)KelvinVoigtTest(t,x,params), tData, ...
        0, options);


    figure(1);
    clf;
    plot(tData, 10^6*xData, 'Marker', 'o', 'LineStyle', 'none');
    hold on;
    plot(tOut, 10^6*xOut);
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
end


