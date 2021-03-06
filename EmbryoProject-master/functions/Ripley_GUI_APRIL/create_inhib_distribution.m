function out = create_inhib_distribution(lambda, delta, x_side, y_side, z_side)
% create_inhib_distribution(lambda, delta, x_side, y_side, z_side)
% Create a distribution with simple inhibition
% Input:
% lambda    - intensity
% side      - size of generated distribution (side of cube)
% delta     - area of inhibition
% Output:
% [X Y Z]
%
% License: RipleyGUI is distributed free under the conditions that
% (1) it shall not be incorporated in software that is subsequently sold; 
% (2) the authorship of the software shall be acknowledged in any publication that uses results generated by the software; 
% (3) this notice shall remain in place in each source file. 

req_pack_intensity = lambda * delta^3 *pi/6;
if(req_pack_intensity > pi/48) % not enough space for distribution
    errordlg('There is not space for so many so sparse events. Try lower intensity or less inhibition','Input error');
    out = [];
    return
end

events = floor(lambda*x_side*y_side*z_side);
pos = [];

while(size(pos,1) < events)
    point = [rand(1)*x_side rand(1)*y_side rand(1)*z_side];
    discard = 0;
    for i = 1:size(pos,1)
        distance = w_distance(pos(i,:),point);
        if distance < delta
            discard = 1; break;
        end
    end
    if discard == 0
        pos = [pos; point];
    end
end

out = pos;