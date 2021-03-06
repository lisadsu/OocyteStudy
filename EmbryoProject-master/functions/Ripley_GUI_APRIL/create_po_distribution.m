function out = create_po_distribution(lambda, x_side, y_side, z_side)
% Create a completely random distribution
% create_po_distribution(lambda, x_side, y_side, z_side)
% lambda  - intensity
% x_side, y_side, z_side    - size of new data set
% Output:
% [X Y Z] positions
%
% License: RipleyGUI is distributed free under the conditions that
% (1) it shall not be incorporated in software that is subsequently sold; 
% (2) the authorship of the software shall be acknowledged in any publication that uses results generated by the software; 
% (3) this notice shall remain in place in each source file. 

events = floor(lambda * x_side * y_side * z_side);
pos = rand(events,3);
pos(:,1) = pos(:,1) * x_side;
pos(:,2) = pos(:,2) * y_side;
pos(:,3) = pos(:,3) * z_side;
out = pos;