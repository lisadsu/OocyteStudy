% Find parameters that solve Kelvin-Voigt Model ODE (KelvinVoigtTest.m)
% x(0) = 0
% 

function f = KelvinVoigtTestFit(tData, xData, F, paramsIn)

k = paramsIn(1);
n = paramsIn(2);

params = [F, k, n];
options = odeset();

% x(0) = F/k;
[tOut, xOut] = ode45( @(t,x)KelvinVoigtTest(t,x,params), tData, ...
    0, options);
% 
% % measure mse
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
% weighting(1) = weighting(1)*10;
f = sum(weighting.*(xOut' - xData).^2);

if k < .01 || n < 0
    
    f = f*2;
    
end