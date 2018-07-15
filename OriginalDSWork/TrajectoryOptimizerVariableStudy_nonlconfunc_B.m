function [c,ceq] = TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type)
if strcmp(fit_type,'exponential')
    num_terms = length(b)/2;
    switch num_terms
        case 1
            r = b(1).*exp(-b(2).*x);
%             ceq = [r(end)-rf;...
%                 b(2)-deg2rad(5)];
        case 2
            r = b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x);
%             ceq = [r(end)-rf;...
%                 b(2)-deg2rad(5);...
%                 b(4)-deg2rad(5)];
        case 3
            r = b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x) + b(5).*exp(-b(6).*x);
%             ceq = [r(end)-rf;...
%                 b(2)-deg2rad(5);...
%                 b(4)-deg2rad(5);...
%                 b(6)-deg2rad(5)];
        case 4
            r = b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x) + b(5).*exp(-b(6).*x) + b(7).*exp(-b(8).*x);
%             ceq = [r(end)-rf;...
%                 b(2)-deg2rad(5);...
%                 b(4)-deg2rad(5);...
%                 b(6)-deg2rad(5);...
%                 b(8)-deg2rad(5)];
        case 5
            r = b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x) + b(5).*exp(-b(6).*x) + b(7).*exp(-b(8).*x) + b(9).*exp(-b(10).*x);
%             ceq = [r(end)-rf;...
%                 b(2)-deg2rad(5);...
%                 b(4)-deg2rad(5);...
%                 b(6)-deg2rad(5);...
%                 b(8)-deg2rad(5);...
%                 b(10)-deg2rad(5)];
        case 6
            r = b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x) + b(5).*exp(-b(6).*x) + b(7).*exp(-b(8).*x) + b(9).*exp(-b(10).*x) + b(11).*exp(-b(12).*x);
%             ceq = [r(end)-rf;...
%                 b(2)-deg2rad(5);...
%                 b(4)-deg2rad(5);...
%                 b(6)-deg2rad(5);...
%                 b(8)-deg2rad(5);...
%                 b(10)-deg2rad(5);...
%                 b(12)-deg2rad(5)];
        otherwise
            disp('Need to create new exponential function in switch statement');
    end
    c = [];
%         ceq = [(r(end)-r(end-1))/x(end) + 0.05;...
%             r(end)-rf];
    ceq = r(end)-rf;
    
elseif strcmp(fit_type,'polynomial')
    num_terms = length(b);
    switch num_terms
        case 1
            r = b(1).*x.^0;
        case 2
            r =b(2).*x.^0 + b(1).*x.^1;
        case 3
            r = b(3).*x.^0 + b(2).*x.^1 + b(1).*x.^2;
        case 4
            r = b(4).*x.^0 + b(3).*x.^1 + b(2).*x.^2 + b(1).*x.^3;
        case 5
            r = b(5).*x.^0 + b(4).*x.^1 + b(3).*x.^2 + b(2).*x.^3 + b(1).*x.^4;
        case 6
            r = b(6).*x.^0 + b(5).*x.^1 + b(4).*x.^2 + b(3).*x.^3 + b(2).*x.^4 + b(1).*x.^5;
        case 7
            r = b(7).*x.^0 + b(6).*x.^1 + b(5).*x.^2 + b(4).*x.^3 + b(3).*x.^4 + b(1).*x.^5 + b(1).*x.^6;
        case 8
            r = b(8).*x.^0 + b(7).*x.^1 + b(6).*x.^2 + b(5).*x.^3 + b(4).*x.^4 + b(3).*x.^5 + b(2).*x.^6 + b(1).*x.^7;
        case 9
            r = b(9).*x.^0 + b(8).*x.^1 + b(7).*x.^2 + b(6).*x.^3 + b(5).*x.^4 + b(4).*x.^5 + b(3).*x.^6 + b(2).*x.^7 + b(1).*x.^8;
        case 10
            r = b(10).*x.^0 + b(9).*x.^1 + b(8).*x.^2 + b(7).*x.^3 + b(6).*x.^4 + b(5).*x.^5 + b(4).*x.^6 + b(3).*x.^7 + b(2).*x.^8 + b(1).*x.^9;
        case 11
            r = b(11).*x.^0 + b(10).*x.^1 + b(9).*x.^2 + b(8).*x.^3 + b(7).*x.^4 + b(6).*x.^5 + b(5).*x.^6 + b(4).*x.^7 + b(3).*x.^8 + b(2).*x.^9 + b(1).*x.^10;
        case 12
            r = b(12).*x.^0 + b(11).*x.^1 + b(10).*x.^2 + b(9).*x.^3 + b(8).*x.^4 + b(7).*x.^5 + b(6).*x.^6 + b(5).*x.^7 + b(4).*x.^8 + b(3).*x.^9 + b(2).*x.^10 + b(1).*x.^11;
    end
    c = [];
    ceq = [r(end)-rf;...
        r(1) - r0];
end
