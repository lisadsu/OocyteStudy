% solve for system parameters using fminsearch
% 
% startParams = [k0 k1 n]

function [paramsOut, fval] = SLSOptimizeParams(tData, xData, startParams, F, plotInput)


[paramsOut, fval, exitflag] = fminsearch(@(paramsIn)SLSTestFit(tData, ...
    xData, F, paramsIn), startParams, ...
    optimset('TolFun', 10^(-14), 'TolX', 10^(-4)));

if (plotInput)
    
    k0 = paramsOut(1);
    k1 = paramsOut(2);
    n = paramsOut(3);

    params = [F, k0, k1, n];
    options = odeset();

    % x(0) = 0;
    [tOut, xOut] = ode45( @(t,x)SLSTest(t,x,params), tData, ...
        F/(k0 + k1), options);

    figure(1);
    clf;
    plot(tData, 10^6*xData, 'Marker', 'o', 'LineStyle', 'none');
    hold on;
    plot(tOut, 10^6*xOut);
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
end


