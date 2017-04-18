% Set error bar width 
% Inputs: errorHandle: handle to errorbar corresponding to bar
%         xCenter: x-coordinate that bar is centered on
%         width: desired errorbar width

function [] = setErrorBar(errorHandle, xCenter, width)

% set error bar width to be width
h = get(errorHandle, 'children');
x = get(h, 'xdata');
x(2) = {[xCenter xCenter NaN xCenter-width xCenter+width ...
    NaN xCenter-width xCenter+width NaN]};
set(h(2), 'xdata', x{2});
set(errorHandle, 'children', h);



