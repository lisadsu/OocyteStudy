% BW = output image from canny edge detector, x and y are bounding box
% around the top corner of the pipette

function [xout yout] = returnCorner(BW, coords)

coords = round(coords);
xmin = coords(1);
xmax = coords(2);
ymin = coords(3);
ymax = coords(4);

BWRightEdge = (BW(ymin:ymax,xmax))';
topRightEdge = find(BWRightEdge == 1, 1);

% yc is increasing going up, so have to invert it
% this is the coordinate of the rightmost end of edge
% xc and yc are referenced to bounding rectangle, not all of BW
xc(1) = xmax - xmin + 1;
yc(1) = ymax - ymin + 1 - topRightEdge + 1;

keepSearching = 1;
currx = xmin + xc(1) - 1;
curry = ymin + topRightEdge - 1;
BWbox = BW(curry-1:curry+1, currx-1:currx+1);

% search directly below first, then go CW around
xsearch = [0 -1 -1 1];% 0 1 1 1];
ysearch = [1 1 0 1];% -1 -1 0 1];

while keepSearching
    
    if (currx - 1 < xmin) || (curry + 1 > ymax)
        keepSearching = 0;
    end
    
    % start below, then go around clockwise
    % once find the next 1, set currx and curry to that, and keep searching
    foundNext = 0;
    i = 1;
    
    while ~foundNext

        if BW(curry + ysearch(i), currx + xsearch(i))
            foundNext = 1;
            currx = currx + xsearch(i);
            xc = [xc (currx - xmin + 1)];
            curry = curry + ysearch(i);
            yc = [yc (ymax - curry + 1)];
        elseif i == 4
            str = 'Couldnt find next';
            keepSearching = 0;
            foundNext = 1;
            msgbox('Problem with contour selected', 'Button', 'warn');
            pause;
        end
        i = i+1;
    end
    
end

% now we have a list, currx, curry of points along the edge
% now find the point where the derivative is closest to 1
% that point is (xout, yout)
pList = [];

for j = 1:length(xc)-6
    
    xSubset = xc(j:j+6);
    ySubset = yc(j:j+6);
    
    if range(xSubset) < 1
        pList = [pList 1000];
    else
        p = polyfit(xSubset, ySubset, 1);
        pList = [pList p(1)];
    end
    
end

pListSub = abs(pList - 1);
bestCoord = find(pListSub == min(pListSub), 1);

bestX = xc(bestCoord);
bestY = yc(bestCoord);

xout = xmin + bestX - 1;
yout = ymax - bestY + 1;
            
            
            
            
            
            
            
            
