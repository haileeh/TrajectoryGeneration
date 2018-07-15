function acrossb = fast_cross(a,b)
% Computes the cross product without the extra MATLAB overhead
%
% Inputs:
% - a:          Vector (3x1)
% - b:          Vector (3x1)
%
% Output:
% - acrossb:    Cross product of a and b
%
% Author: Christopher M. Pong
%
%#codegen

acrossb = [
    a(2)*b(3) - a(3)*b(2)
    a(3)*b(1) - a(1)*b(3)
    a(1)*b(2) - a(2)*b(1)];

end
