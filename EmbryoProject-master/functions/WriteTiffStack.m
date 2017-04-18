
% write 3D double matrix as 16-bit TIFF stack

function [] = WriteTiffStack(matrix, filename)

for i = 1:size(matrix,3)
    
    if (i == 1)
        imwrite(uint16(floor((2^16)*matrix(:,:,i))), filename, 'TIFF');
    else
        imwrite(uint16(floor((2^16)*matrix(:,:,i))), filename, 'TIFF', ...
            'writemode', 'append');
    end
    
end


