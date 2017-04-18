% t and aspiration depth should start at the data point right before the
% cell enters the pipette (thus the first data point should be at 0, and
% the rest should be higher)

% aspiration depth in KelvinFit2 needs to be in meters, not pixels

function [xfine yfit k0 k1 n0 n1_inv F0 tau fval] = ...
    KelvinFit3(t, aspiration_depth, Fin, plot_input, startParams, colorIn)

if nargin < 6
    colorIn = [0 0 1];
end

F0 = Fin; % approx pressure*area
start_params = startParams;
% start_params(1) = .02; % k0
% start_params(2) = .067; % k1
% start_params(3) = .2; % tau
% start_params(4) = .5; % n1_inv (slope of creep)

% add a 1/(k0 + k1) offset to represent instant elongation, start the
% exponential fit from there.
xdata = t(2:end);
% ydata = [F0/(start_params(1) + start_params(2)) aspiration_depth(2:end)]0;
ydata = aspiration_depth(1:end-1);

options = optimset('TolFun', 10^-14);
[best_params fval exitflag] = fminsearch(@(params)KelvinModel2(params, xdata, ...
    ydata, F0), start_params, options);

k0 = best_params(1);
k1 = best_params(2);
tau = best_params(3);
n1_inv = best_params(4);

n0 = tau*(k0*k1)/(k0 + k1);

% recalculate fit and plot, with finer sampling
xfine = linspace(min(xdata), max(xdata), 1000);
yfit = F0/k1*(1 - k0/(k0+k1)*exp(-xfine/tau)) + xfine*F0*n1_inv;

% plot, first point of aspiration depth is replaced by the height after
% instant elongation

if (plot_input == 1)
    % figure;
    plot(t, 10^6*[F0/(k0 + k1) aspiration_depth(1:end-1)],'ob', 'Color', colorIn);
    hold on;
    % this needs to be shifted to account for truncation during fitting
    plot(xfine, 10^6*yfit, 'Color', colorIn);
    set(gca, 'FontSize', 14);
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
    title('Aspiration Depth into Micropipette');
    % hold off;
end
