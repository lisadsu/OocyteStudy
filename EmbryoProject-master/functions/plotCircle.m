function plotCircle(p1, p2, p3)


p1 = [448 366]
p2 = [790 418]
p3 = [1136 326]

[c r] = calcCircle(p1, p2, p3)

[xC, yC] = circlepoints(r);
xC = (xC+c(1));
yC = (yC+c(2));

outOfBounds = (xC < (startInd*scaleFactor)) | (xC > (endInd*scaleFactor))...
    | (yC < 1);% | (xC < 1) | ...
%    (xC > size(RF,2));

xC(outOfBounds) = NaN;
yC(outOfBounds) = NaN;

xC = xC(~isnan(xC));
yC = yC(~isnan(yC));

[xS, I] = sort(xC, 'ascend');
yS = yC(I);

hold on;
plot(yS, xS, 'b-', 'Color', [1 1 0], 'LineWidth', 3);


end



