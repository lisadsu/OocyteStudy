% Find parameters that solve MOdified SLS Model ODE (SLSTest.m)
% x(0) = F / (k1 + k0)
% 

function f = FiveParamTestFit(tData, xData, F, paramsIn)

k0 = paramsIn(1);
k1 = paramsIn(2);
k2 = paramsIn(3);
n0 = paramsIn(4);
n1 = paramsIn(5);

params = [F, k0, k1, k2, n0, n1];
options = odeset();

% x(0) = F/(k0 + k1);
x0 = F/(k0 + k1 + k2);

[tOut, xOut] = ode45( @(t,x)FiveParamTest(t,x,params), tData, ...
    [x0 0 x0 0 x0], options);

% f = sum((xOut(:,4)' - xData).^2);

% place extra weight on good fit in first data points
weighting = 1*ones(1, length(xData));

if length(xData) > 20
    weighting(1:10) = 2;
    weighting(1) = 10;
else
    weighting(1:length(xData)-2) = 1;
end

% make the first point extra important that the line go thru
% weighting(1) = weighting(1)*10;
f = sum(weighting.*(xOut(:,5)' - xData).^2);

if k0 < 0.0001 || k1 < 0.0001 || k2 < 0.0001 || n0 < 0.0001 || n1 < 0.0001
    
    f = f*100000;
    
end
