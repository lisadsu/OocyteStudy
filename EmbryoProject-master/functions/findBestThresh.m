function threshOut = findBestThresh(ROIframes, pCenter, pRange, wMaxInd, wMaxH)

% find best canny threshold to just barely not detect anything inside empty
% pipette

threshList = linspace(1,0,50);
threshList = threshList(2:end-1);
sumList = zeros(1,length(threshList));

for i = 1:length(threshList)
    c1 = edge(ROIframes(:,:,1), 'canny', threshList(i));
    c1 = c1(pCenter-pRange:pCenter+pRange,3*wMaxInd:wMaxH);
    sumList(i) = sum(sum(c1));
end

figure, plot(diff(sumList))

bestThresh = find(diff(sumList) < 10, 1, 'last');
% bestThresh = 45
threshOut = threshList(bestThresh)

cThresh = edge(ROIframes(:,:,1), 'canny', threshOut);
cThresh = cThresh(pCenter-pRange:pCenter+pRange,3*wMaxInd:wMaxH);

cThresh2 = edge(ROIframes(:,:,1), 'canny', threshList(bestThresh + 1));
cThresh2 = cThresh2(pCenter-pRange:pCenter+pRange,3*wMaxInd:wMaxH);
    
figure, imshow(cThresh);
figure, imshow(cThresh2);

pause(2);
