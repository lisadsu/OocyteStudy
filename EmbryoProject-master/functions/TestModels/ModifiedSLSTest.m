% Modified SLS : spring in parallel with spring + dashpot, plus extra
% dashpot
%
% x'' = k0*k1*(F - n1 * x') / (n0 * n1 * (k0 + k1))
% x = [x1', x2', x3', x']
% x1' = F/n1
% x2' = -k0*k1*x2/n0/(k0 + k1)
% x3' = k0*x2/n0

function dx = ModifiedSLSTest(t,x,params)

F = params(1);
k0 = params(2);
k1 = params(3);
n0 = params(4);
n1 = params(5);

dx = zeros(4,1); 
dx(1) = F / n1;
dx(2) = -1*k0*k1*x(2)/(n0*(k0 + k1));
dx(3) = k0*x(2)/n0;
dx(4) = dx(1) + dx(2) + dx(3);