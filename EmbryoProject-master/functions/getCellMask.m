% Function: getCellMask
% Usage: [mask, mR, mX, mY] = getCellMask(ROI)
%
% ROI should be a 2D image with a single circular cell in it
% This function asks user to either manually make a mask with
%   roipoly, or will use a Hough transform to automatically detect
%   a circle of the appropriate size.
%
%

function [mask, mR, mX, mY] = getCellMask(ROI, radii)

while (true)
    autoMask = input('Automatically detect cell body? (1 for yes, 2 for no). ');
    if (autoMask ~= 1 && autoMask ~= 2)
        disp('Invalid input.');
    else
        break;
    end
end

mX = NaN;
mY = NaN;
mR = NaN;
mask = NaN;

while (isnan(mX))
    
    if autoMask == 1
        % auto mode
        
        clf;
        imshow(imadjust(ROI));
        
        im = imfilter(imadjust(ROI), fspecial('gaussian', [20 20], 10));
        im = medfilt2(im, [10 10]);
        
        % automatic thresholding
        lowThresh = graythresh(im(10:end-10,10:end-10));
        im(im < lowThresh) = 0;
        im(im > lowThresh) = 1;
 
        sliceEdge = edge(im, 'canny');%, [.005 .05]);
        
        % dilate and multiply by local intensity
        edgesDilate = imdilate(sliceEdge, strel('disk', 20));
        edgesBright = edgesDilate .* im;

        h = circle_hough(edgesBright, radii, 'same', 'normalise');
        
        % Find some peaks in the accumulator
        nHoodXY=51;
        nHoodR=51;
        nPeaks=1;
        peaks = circle_houghpeaks(h, radii, 'nhoodxy', nHoodXY, 'nhoodr',...
            nHoodR, 'npeaks', nPeaks, 'smoothxy', 4);
        
        if ~isempty(peaks)
            [x, y] = circlepoints(peaks(3));
            x = x+peaks(1);
            y = y+peaks(2);
            
            hold on;
            plot(x, y, 'g-', 'LineWidth', 5);
            
            mX = peaks(1); % mX pixels to the right
            mY = peaks(2); % mY pixels down
            mR = peaks(3);
            
            % make mask
            [iC, jC] = meshgrid(1:size(ROI,2), 1:size(ROI,1));
            D = sqrt((iC - mX).^2 + (jC - mY).^2);
            
            mask = zeros(size(ROI,1), size(ROI,2));
            mask(D < mR) = 1;
            mask = (mask == 1);
                        
            while (true)
                exitProgram = input('Satisfied with cell detection? (1 for yes, 2 for no). ');
                if (exitProgram ~= 1 && exitProgram ~= 2)
                    disp('Invalid input.');
                elseif (exitProgram == 1)
                    return;
                else
                    autoMask = 2; % if user doesn't like mask, switch to manual mode
                    mX = NaN;
                    mY = NaN;
                    mR = NaN;
                    mask = NaN;
                    break;
                end
            end
            
        else
            % if no peak was found, switch to manual mode
            autoMask = 2;
        end
        
    else
        % manual mode
        clf;
        imshow(imadjust(ROI));
        disp('select cell outline with polygon');
        title('select cell outline with polygon');
        BW = roipoly;
        
        % estimate x, y, r of region
        stats = regionprops(BW, 'Centroid', 'BoundingBox');
        mX = stats.Centroid(1);
        mY = stats.Centroid(2);
        mR = mean(stats.BoundingBox(3:4))/2;
        
        [x, y] = circlepoints(mR);
        x = x+mX;
        y = y+mY;
        
        hold on;
        plot(x, y, 'g-', 'LineWidth', 5);
        
        % make mask
        [iC, jC] = meshgrid(1:size(ROI,2), 1:size(ROI,1));
        D = sqrt((iC - mX).^2 + (jC - mY).^2);
        
        mask = zeros(size(ROI,1), size(ROI,2));
        mask(D < mR) = 1;
        mask = (mask == 1);
        
        while (true)
            exitProgram = input('Satisfied with cell outline? (1 for yes, 2 for no). ');
            if (exitProgram ~= 1 && exitProgram ~= 2)
                disp('Invalid input.');
            elseif (exitProgram == 1)
                return;
            else
                autoMask = 2; % if user doesn't like mask, switch to manual mode
                mX = NaN;
                mY = NaN;
                mR = NaN;
                mask = NaN;
                break;
            end
        end
        
    end
    
end









