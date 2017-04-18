% Kelvin-Voigt body : spring + dashpot in parallel
%
% x' = F/n

function dx = KelvinVoigtTest(t,x,params)

F0 = params(1);
k = params(2);
n = params(3);

dx = F0/n - k*x/n;