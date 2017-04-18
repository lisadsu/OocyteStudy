% Xref is just 0:.01:1
% Yi is Y at the values in Xref instead of X
% can't use interp1 because X and Y have duplicate values

function Yi = interpForRoc(Xref, X, Y)


% get rid of duplicate values - just pick average Y value
Xdiff = diff(X);
currInd = 1;
Xnew = [];
Ynew = [];

while currInd < length(Xdiff) - 1
    
    if Xdiff(currInd) > 0
        
        % if not in a vertical region, keep current X and Y values
        Xnew = [Xnew X(currInd)];
        Ynew = [Ynew Y(currInd)];
        currInd = currInd + 1;
        
    else
        
        % keep just one X value from the vertical region
        % keep the average Y value from the vertical region
        XdiffRest = Xdiff(currInd:end);
        vertEnd = find(XdiffRest > 0, 1);
        
        Xnew = [Xnew X(currInd)];
        if ~isempty(vertEnd) && currInd + vertEnd < length(Y)
            Ynew = [Ynew mean(Y(currInd:currInd+vertEnd))];
        else
            Ynew = [Ynew mean(Y(currInd:end))];
        end
        
        currInd = currInd + vertEnd;
        
    end
        
end

if min(length(Xnew), length(Ynew)) > 2
    Yi = interp1(Xnew, Ynew, Xref, 'linear', 'extrap');
else
    Yi = [0 ones(1,length(Xref)-1)];
end