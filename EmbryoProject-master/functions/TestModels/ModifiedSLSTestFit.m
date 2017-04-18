% Find parameters that solve MOdified SLS Model ODE (SLSTest.m)
% x(0) = F / (k1 + k0)
% 

function f = ModifiedSLSTestFit(tData, xData, F, paramsIn)

k0 = paramsIn(1);
k1 = paramsIn(2);
n0 = paramsIn(3);
n1 = paramsIn(4);

params = [F, k0, k1, n0, n1];
options = odeset();

% x(0) = F/(k0 + k1);
% x'(0) = 
x0 = F/(k0 + k1);
dx0 = F/n1 + F/(k0 + k1) + F*(k0^2)/(n0 * (k0 + k1)^2);

[tOut, xOut] = ode45( @(t,x)ModifiedSLSTest(t,x,params), tData, ...
    [0 x0 0 x0], options);

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
f = sum(weighting.*(xOut(:,4)' - xData).^2);

if k0 < 0 || k1 < 0 || n0 < 0 || n1 < 0
    
    f = f*2;
    
end
