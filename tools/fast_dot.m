function adotb = fast_dot(a,b)
% Computes the dot product without the extra MATLAB overhead
%
% Inputs:
% - a:      Vector
% - b:      Vector
%
% Output:
% - adotb:  Dot product of a and b
%
% Author: Christopher M. Pong
%
%#codegen

adotb = a.'*b;

end
