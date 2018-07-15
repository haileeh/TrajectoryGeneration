function q_A_B = conj_quat(q_B_A)
% Calculates the inverse/conjugate of the quaternion
%
% Input:
% - q_B_A:  Transformation to frame B from A (4x1)
%
% Output:
% - q_A_B:  Transformation to frame A from B (4x1)
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

q_A_B = [-q_B_A(1:3); q_B_A(4)];

end
