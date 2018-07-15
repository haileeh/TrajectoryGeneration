function dcm = dcm_from_quat(quat)

q1 = quat(1); q2 = quat(2); q3 = quat(3); q0 = quat(4);

dcm = [(q0^2+q1^2-q2^2-q3^2)    2*(q1*q2+q0*q3) 2*(q1*q3-q0*q2);...
        2*(q1*q2 - q0*q3)       (q0^2-q1^2+q2^2-q3^2)   2*(q2*q3+q0*q1);...
        2*(q1*q3+q0*q2)         2*(q2*q3-q0*q1)     (q0^2-q1^2-q2^2+q3^2)];

end