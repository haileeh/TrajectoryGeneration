function quat3 = mult_quat(quat1,quat2)
% For use, when scalar is 4th element

r1 = quat1(1); r2 = quat1(2); r3 = quat1(3); r0 = quat1(4);
p1 = quat2(1); p2 = quat2(2); p3 = quat2(3); p0 = quat2(4);

n1 = r0*p1 + r1*p0 -r2*p3 + r3*p2;
n2 = r0*p2 + r1*p3 + r2*p0 - r3*p1;
n3 = r0*p3 - r1*p2 + r2*p1 + r3*p0;
n0 = -r1*p1 - r2*p2 - r3*p3 + r0*p0;

quat3 = [n1 n2 n3 n0]';

end