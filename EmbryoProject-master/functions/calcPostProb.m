function [xfit, yfit] = calcPostProb(xdata, ydata)

% calculate fit for adjViability vs log(distList) just by binning
% include smoothing since some bins might have 0 observations

xfit = linspace(min(xdata), max(xdata), 500);
yfit = zeros(1,length(xfit));

for i = 1:length(xfit) - 1
    
    currMin = xfit(i);
    currMax = xfit(i+1);
    
    if isempty(ydata((xdata <= currMax) & (xdata >= currMin)))
        % if there's no data in this interval, set y to NaN
        % or maybe expand interval instead?
        yfit(i) = NaN;
    elseif isempty(ydata((xdata <= currMax) & (xdata >= currMin) & (ydata == 1)))
        % if there are no 1s, then set yfit to 0
        yfit(i) = 0;
    elseif isempty(ydata((xdata <= currMax) & (xdata >= currMin) & (ydata == 0)))
        % if there are no 0s, then set yfit to 1
        yfit(i) = 1;
    else
        yfit(i) = mean(ydata(xdata < currMax & xdata > currMin));
    end

end

% get rid of NaNs

for i = 1:length(xfit) - 1
    
   if isnan(yfit(i))
       
      % find first non-nan before and after current bin and average
      
      lowNonNan = find(~isnan(yfit(1:i-1)), 1, 'last');
      highNonNan = find(~isnan(yfit(i+1:end)), 1, 'first');
      
      yfit(i) = mean([yfit([lowNonNan highNonNan])]);
       
   end
    
end


