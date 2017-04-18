% convert video frames to TIFF series
% input is the path to the .avi file you want to convert (with backslashes,
%   not forward slashes)
% you can also call the function with no argument and a dialog box will pop
%   up so you can select a file

function [] = convertVideoToImageSeries(input)

if nargin < 1
    [filename, path] = uigetfile({'*.avi'}, 'Video File Select');
    input = [path filename];
end

% read in video
obj = VideoReader(input);
frames = squeeze(read(obj));

% make output directory
[pathstr, name, ~] = fileparts(input);
outputDir = [pathstr, '\', name, '_series\'];

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for i = 1:obj.NumberOfFrames
    
    % write out TIFF series with number field 3 chars wide
    % accept both grayscale and color inputs
    if length(size(frames)) == 3
        imwrite(frames(:,:,i), [outputDir name '_' num2str(i, '%0.3d') '.tif']);
    else
        imwrite(frames(:,:,:,i), [outputDir name '_' num2str(i, '%0.3d') '.tif']);
    end
    
end












