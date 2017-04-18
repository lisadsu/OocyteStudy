function f = KelvinModel2(params, xdata, ydata, F0)

k0 = params(1);
k1 = params(2);
tau = params(3);
n1_inv = params(4);

ycalc = F0/k1*(1 - k0/(k0+k1)*exp(-xdata/tau)) + xdata*F0*n1_inv;

% adjust height of first point of ydata depending on what the current
% values of k0 and k1 are
% ydata(1) = F0/(k0 + k1);

% the fit in the curved part is more important than the fit at the straight
% part, which tends to dominate since it's most of the data
% weighting = exp(-xdata/1);
weighting = .1*ones(1, length(ydata));

if length(ydata) > 10
    weighting(1:10) = 1;
else
    weighting(1:length(ydata)-2) = 1;
end

% make the first point extra important that the line go thru
weighting(1) = weighting(1)*10;

% f = sum(weighting .* ((ycalc - ydata).^2));
f = sum(weighting.*(ycalc - ydata).^2);

if k0 < 0 || k1 < 0 || n1_inv < 0 || tau < 0
    
    f = f*10;
    
end

