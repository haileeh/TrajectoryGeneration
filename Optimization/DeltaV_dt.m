function DeltaV = DeltaV_dt(dt_list,omega,r0,rf,inertia_ratio,exp_coefficients)

if nargin < 6
    work_with_exponentials = 0;
else
    work_with_exponentials = 1;
end

%Clear the data vectors
r_SVC_INT = [r0; 0; 0];  % m
r_SVC_TAR = r_SVC_INT;  % ASSUMPTION
r0 = norm(r_SVC_INT);

I_TAR = diag(inertia_ratio);
omega_TAR_INT = deg2rad(omega); %rate of target in inertial frame

v_SVC_INT = fast_cross(omega_TAR_INT,r_SVC_INT); %linear velocity of servicer relative to static target in inertial frame
% Make sure this is consistent with other deltaV function
q_TAR_INT = [0; 0; 0; 1];   %quaternion from inertial frame to target frame
q_UDP_TAR = [0; 0; 0; 1];   %quaternion from target frame to inertial frame
% q_UDP_TAR = [0.0019;   -0.0436;    0.0436;    0.9981];
% q_UDP_TAR = [0 0.7071 0.7071 0]';

% Initialize variables
initialize_zeros 
%rdot_SVC_TAR = [0.2 0 0]'; %[0.01 0 0]';
num_iter = length(dt_list);

if work_with_exponentials
    rad_dist = exp_coefficients(1)*exp(-exp_coefficients(2)*t) + ...
       exp_coefficients(3)*exp(-exp_coefficients(4)*t);  
else % i.e. for full optimization
    rad_dist = logspace(log10(r0),log10(rf),num_iter+1);
end
r = r0;
for i = 1:num_iter
    if r(i) < rf
        break;
    end
    
    dt = dt_list(i);
    tumbling_dynamics
    
    if work_with_exponentials == 0
        %
    else
        %
    end
    
    % Calculate accelerations
    acc_DV_calculations
    
    % Normalize r to compare against rf
    r(i+1) = norm(r_SVC_TAR(:,i+1));

    t(i+1) = t(i) + dt;
end % end loop of simulating time until end of spiral
% DV_total = smooth(t,DV_total,0.1,'rloess');
assignin('base','funcT',t)
assignin('base','funcR',rad_dist);

DV_TOT = DV_total(end); 

DeltaV = DV_TOT;

end