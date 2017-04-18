% solve for system parameters using fminsearch
% 
% startParams = [k0 k1 n]

function [paramsOut, fval] = FiveParamOptimizeParams(tData, xData, startParams, F, plotInput)


[paramsOut, fval, exitflag] = fminsearch(@(paramsIn)FiveParamTestFit(tData, ...
    xData, F, paramsIn), startParams, ...
    optimset('TolFun', 10^(-16), 'TolX', 10^(-6), 'MaxFunEvals', 100000));


if (plotInput)
    
    k0 = paramsOut(1);
    k1 = paramsOut(2);
    k2 = paramsOut(3);
    n0 = paramsOut(4);
    n1 = paramsOut(5);

    params = [F, k0, k1, k2, n0, n1];
    paramsOut
    fval
    options = odeset();

    % x(0) = F/(k0 + k1);
    x0 = F/(k0 + k1 + k2);

    [tOut, xOut] = ode45( @(t,x)FiveParamTest(t,x,params), tData, ...
        [x0 0 x0 0 x0], options);


    figure(1);
    clf;
    plot(tData, 10^6*xData, 'Marker', 'o', 'LineStyle', 'none');
    hold on;
    plot(tOut, 10^6*xOut(:,5));
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
end


