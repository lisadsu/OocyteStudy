% solve for system parameters using fminsearch
% 
% startParams = [k0 k1 n]

function [paramsOut, fval] = ModifiedSLSOptimizeParams(tData, xData, startParams, F, plotInput)


[paramsOut, fval, exitflag] = fminsearch(@(paramsIn)ModifiedSLSTestFit(tData, ...
    xData, F, paramsIn), startParams, ...
    optimset('TolFun', 10^(-16), 'TolX', 10^(-6), 'MaxFunEvals', 1000));


if (plotInput)
    
    k0 = paramsOut(1);
    k1 = paramsOut(2);
    n0 = paramsOut(3);
    n1 = paramsOut(4);

    params = [F, k0, k1, n0, n1];
    options = odeset();

    % x(0) = F/(k0 + k1);
    % x'(0) = ;
    x0 = F/(k0 + k1);
    dx0 = F/n1 + F/(k0 + k1) + F*(k0^2)/(n0 * (k0 + k1)^2);

    [tOut, xOut] = ode45( @(t,x)ModifiedSLSTest(t,x,params), tData, ...
        [0 x0 0 x0], options);

    figure(1);
    clf;
    plot(tData, 10^6*xData, 'Marker', 'o', 'LineStyle', 'none');
    hold on;
    plot(tOut, 10^6*xOut(:,4));
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
    
end


