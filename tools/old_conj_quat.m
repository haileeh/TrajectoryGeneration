function c_quat = conj_quat(quat)
% For use with quaternion with scalar as last element

r1 = quat(1); r2 = quat(2); r3 = quat(3); r4 = quat(4);

c_quat = [-r1 -r2 -r3 r4]';

end