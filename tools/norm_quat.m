function q = norm_quat(q)
% Normalizes and properizes the quaternion
%
% Input:
% - q:  Quaternion to be normalized and properized (4x1)
%
% Output:
% - q:  Normalized and properized quaternion (4x1)
%
% Note:
% - Assumes the quaternion has the vector first and scalar last
%
% Author: Christopher M. Pong
%
%#codegen

% Normalize the quaternion
q = q/norm(q);

% Properize the quaternion (ensure that rotation is less than 180 degrees)
if q(4) < 0
    q = -q;
end

end
