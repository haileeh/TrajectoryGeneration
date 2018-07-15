function q_C_A = mult_quat(q_C_B,q_B_A)
% Multiplies two quaternions
%
% Inputs:
% - q_C_B:  Transformation to frame C from B (4x1)
% - q_B_A:  Transformation to frame B from A (4x1)
%
% Output:
% - q_C_A:  Transformation to frame C from A (4x1)
%
% Reference:
% - Breckenridge, Quaternions-Proposed Standard Conventions, 31 Oct. 1979,
%   JPL Interoffice Memorandum.
%
% Note:
% - Assumes the quaternion has the vector first and scalar last
%
% Author: Christopher M. Pong
%
%#codegen

% Extract scalar part of the quaternions
s_C_B = q_C_B(4);
s_B_A = q_B_A(4);

% Extract the vector part of the quaternions
v_C_B = q_C_B(1:3);
v_B_A = q_B_A(1:3);

% Perform the multiplication
q_C_A = [
    s_C_B*v_B_A + s_B_A*v_C_B - fast_cross(v_C_B,v_B_A)
    s_C_B*s_B_A - fast_dot(v_C_B,v_B_A)];

end
