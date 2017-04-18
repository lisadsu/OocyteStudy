% solve for system parameters using fminsearch
% 
% startParams = [n k]

function [paramsOut, fval] = MaxwellOptimizeParams(tData, xData, startParams, F, plotInput)


[paramsOut, fval, exitflag] = fminsearch(@(paramsIn)MaxwellTestFit(tData, ...
    xData, F, paramsIn), startParams, ...
    optimset('TolFun', 10^(-9), 'TolX', 10^(-4), 'MaxFunEvals', 1000));

if (plotInput)

    k = paramsOut(1);
    n = paramsOut(2);

    params = [F, n];
    options = odeset();

    % x(0) = F/k;
    [tOut, xOut] = ode45( @(t,x)MaxwellTest(t,x,params), tData, ...
        F/k, options);


    figure(1);
    clf;
    plot(tData, 10^6*xData, 'Marker', 'o', 'LineStyle', 'none');
    hold on;
    plot(tOut, 10^6*xOut);
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
end



