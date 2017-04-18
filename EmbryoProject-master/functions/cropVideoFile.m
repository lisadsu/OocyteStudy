% read in video file
% select subROI
% and save that as its own video file

clear all;
close all;

movpath = 'C:\Users\Livia\Desktop\Human embryo map\Embryo movies\Processed\';
startDir = pwd;
cd(movpath);
filename = uigetfile('*.avi');
cd(startDir);

% read in video, keep only first color channel
obj = VideoReader([movpath filename]);
frameRate = obj.FrameRate;
frames = read(obj);
grayframes = squeeze(frames(:,:,1,:));
% clear frames;

%% save multiple sub-ROIs from each loaded video file

outputname = 'BE2_D3';

% display a single frame for the user to crop
figure, imshow(grayframes(:,:,1));
title('select ROI to crop');
coord = round(getrect);
coords = round([coord(1) coord(1)+coord(3) coord(2) coord(2)+coord(4)]);
ROIframes = grayframes(coords(3):coords(4), coords(1):coords(2), :);
% clear grayframes;

objOut = VideoWriter([movpath outputname '.avi'], 'Grayscale AVI');
open(objOut);

for i = 1:size(ROIframes,3)
    writeVideo(objOut, ROIframes(:,:,i));
end

close(objOut);
close all;


