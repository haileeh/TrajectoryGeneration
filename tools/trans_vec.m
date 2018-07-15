function vec_B = trans_vec(q_B_A,vec_A)
% Transforms a vector using a quaternion from frame A to B
%
% Inputs:
% - vec_A:      Vector in frame A (3x1)
% - q_B_A:      Transformation to frame B from A (4x1)
%
% Output:
% - vec_B:      Vector in frame B (3x1)
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

% Perform quaternion multiplication
temporary = mult_quat(mult_quat(q_B_A,[vec_A;0]),conj_quat(q_B_A));

% Extract the vector from the temporary quaternion
vec_B = temporary(1:3);

end
