% auto adjust contrast and brightness through stack to 

function imageMatrixGadjust = AdjustContrastThroughStack(imageMatrixG, figHandle, subplotNum, subplotSize)

imageMatrixGadjust = imageMatrixG;

% get image value range for top few and bottom few images in stack, create
% ramp
nSlices = size(imageMatrixG,3);
sliceList = 1:nSlices;
avgSignal = mean(squeeze(mean(imageMatrixG, 1)),1);
nonzeroSignal = (avgSignal > 0);
if (sum(nonzeroSignal) == 0)
    return;
end

% fit line to nonzero signal
p = polyfit(sliceList(nonzeroSignal), avgSignal(nonzeroSignal), 1); 
yfit = polyval(p, sliceList(nonzeroSignal));
yfit(yfit < 0) = 0;

signalRamp = avgSignal;
signalRamp(nonzeroSignal) = yfit;
refMax = max(avgSignal); % reference is brightest point in average signal
figure(figHandle);
subplot(subplotSize,subplotSize,subplotNum);
plot(signalRamp);
hold on;
plot(avgSignal);
ylim([0 .05]);

% now adjust contrast using signal ramp
for i = 1:nSlices
    
    currMax = max(max(imageMatrixG(:,:,i)));    
    avgMax = signalRamp(i);
    
    if currMax > 0 && avgMax > 0
        
        imageMatrixGadjust(:,:,i) = imadjust(imageMatrixG(:,:,i), ...
            [0 avgMax], [0 refMax]);
        
    end
    
end












