% Maxwell body : spring + dashpot in series
%
% x' = F/n

function dx = MaxwellTest(t,x,params)

F0 = params(1);
n = params(2);
dx = F0 / n;