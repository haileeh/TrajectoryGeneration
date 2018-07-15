function [c,ceq] = TrajectoryOptimizerVariableStudy_nonlconfunc_GLF(x,b)
r = b(5) + (b(6)-b(5))./(b(1)+b(3)*exp(-b(4).*x)).^(1/b(2));
c = [];
% ceq = [r(1)-10; r(end)-1];
ceq = [b(5) + (b(6)-b(5))./(b(1)+b(3)*exp(-b(4).*x(1))).^(1/b(2)) - 10;...
    b(5) + (b(6)-b(5))./(b(1)+b(3)*exp(-b(4).*x(end))).^(1/b(2)) - 1];