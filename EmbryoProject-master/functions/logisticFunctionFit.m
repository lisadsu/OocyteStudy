function fval = logisticFunctionFit(xdata, ydata, params)

numViable = length(ydata(ydata == 1));
numNonviable = length(ydata(ydata == 0));

% fit data to A / (1 + B * exp(-C(x-D)))
A = params(1);
B = params(2);
C = params(3);
D = params(4);

yfit = A * (1 + B * (exp(-1*C*(xdata - D)))) .^ (-1);

% apply weighting since there are unequal numbers of viable and nonviable

fval = sum((yfit - ydata).^2);
% diffViable = yfit(ydata == 1) - ydata(ydata == 1);
% diffNonviable = yfit(ydata == 0) - ydata(ydata == 0);
% 
% fval = sum(diffViable .^2) / numViable + sum(diffNonviable.^2) / numNonviable;