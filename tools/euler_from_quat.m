function [yaw, pitch, roll] = euler_from_quat(quat)
% weird order

q1 = quat(1); q2 = quat(2); q3 = quat(3); q0 = quat(4);

sinr = 2*(q0*q1 + q2*q3);
cosr = 1-2*(q1^2+q2^2);
roll = atan2(sinr,cosr);

sinp = 2*(q0*q2-q3*q1);
pitch = asin(sinp);

siny = 2*(q0*q3 + q1*q2);
cosy = 1-2*(q2^2+q3^2);
yaw = atan2(siny,cosy);

end