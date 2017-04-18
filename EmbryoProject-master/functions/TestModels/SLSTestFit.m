% Find parameters that solve SLS Model ODE (SLSTest.m)
% x(0) = F / (k1 + k0)
% 

function f = SLSTestFit(tData, xData, F, paramsIn)

k0 = paramsIn(1);
k1 = paramsIn(2);
n = paramsIn(3);

params = [F, k0, k1, n];
options = odeset();

% x(0) = F/k;
[tOut, xOut] = ode45( @(t,x)SLSTest(t,x,params), tData, ...
    F/(k0 + k1), options);

% f = sum((xOut' - xData).^2);

% place extra weight on good fit in first data points
weighting = 1*ones(1, length(xData));

if length(xData) > 20
    weighting(1:10) = 2;
    weighting(1) = 10;
else
    weighting(1:length(xData)-2) = 1;
end


% make the first point extra important that the line go thru
weighting(1) = weighting(1)*10;
f = sum(weighting.*(xOut' - xData).^2);

if k0 < .01 || k1 < .01 || n < 0
    
    f = f*2;
    
end

