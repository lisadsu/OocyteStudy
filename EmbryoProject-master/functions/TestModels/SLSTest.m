% SLS : spring in parallel with spring + dashpot
%
% x' = k0 / (k0+k1) / n * (F - k1 * x)

function dx = SLSTest(t,x,params)

F0 = params(1);
k0 = params(2);
k1 = params(3);
n = params(4);

dx = k0 * (F0 - k1 * x) / (n * (k0 + k1));