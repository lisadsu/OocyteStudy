% 5-parameter model: spring in parallel with two spring+dashpot series
%
% x = [x1', x2', x3', x4', x']
% x1' = F/n1
% x2' = -k0*k1*x2/n0/(k0 + k1)
% x3' = k0*x2/n0

function dx = FiveParamTest(t,x,params)

F = params(1);
k0 = params(2);
k1 = params(3);
k2 = params(4);
n0 = params(5);
n1 = params(6);

dx = zeros(4,1); 
dx(1) = (k2*n1*(k1+k2)*(F - k1*x(3) - k2*x(5)) - k1*n0*k2*(F - k0*x(1) - k2*x(5))) / ...
    (n0*n1*(k0*k1 - (k0 + k2)*(k1 + k2)));
dx(2) = (F - k1*x(3) - k2*x(5))/n0;
dx(3) = -1*((k0+k2)*dx(1) + k2*dx(2))/k1;
dx(4) = (F - k0*x(1) - k2*x(5))/n1;
dx(5) = dx(1) + dx(2);