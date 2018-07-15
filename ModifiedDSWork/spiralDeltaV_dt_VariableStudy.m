function DeltaV = spiralDeltaV_dt_VariableStudy(dt_list,omega,R0,RF,inertia_ratios,exp_coefficients)

if nargin < 6
    work_with_exponentials = 0;
else
    work_with_exponentials = 1;
    exp_coefficients;
end

OMEGA = size(omega,2); %deg/sec
INERTIA = size(inertia_ratios,2);

for omega_iter = 1:length(OMEGA)
    for r0_iter = 1:length(R0)
        for rf_iter = 1:length(RF)
            for INERTIA_idx = 1:length(INERTIA)
                
                %Set the trajectory parameters
                rf = RF(rf_iter);
                r0 = R0(r0_iter);
                omega_init = norm(omega);
                
                %Clear the data vectors
                r_SVC_INT = [r0; 0; 0];  %m
                r_SVC_TAR = r_SVC_INT;  % ASSUMPTION
                r0 = norm(r_SVC_INT);
               
                I_TAR = diag(inertia_ratios(:,INERTIA_idx));
                omega_TAR_INT = deg2rad(omega); %rate of target in inertial frame
                
                v_SVC_INT = fast_cross(omega_TAR_INT,r_SVC_INT); %linear velocity of servicer relative to static target in inertial frame
                q_TAR_INT = [0; 0; 0; 1];   %quaternion from inertial frame to target frame
                q_UDP_TAR = [0; 0; 0; 1];   %quaternion from target frame to inertial frame
%                 q_UDP_TAR = [0.0019;   -0.0436;    0.0436;    0.9981];
%                 q_UDP_TAR = [0 0.7071 0.7071 0]';
                
                num_iter = length(dt_list); % num_iter = 1000000;
                
                % Initialize variables
                t = 0;
                q_UDP_INT = zeros(4,1);  q_SVC_INT = zeros(4,1);
                M_TAR = zeros(3,1); tau_TAR = zeros(3,1);
                omega_dot_TAR_INT = zeros(3,1);
                q_dot_TAR_INT = zeros(4,1); r_TEMP = zeros(3,1);
                rdot_SVC_TAR = zeros(3,length(dt_list));
                
                if work_with_exponentials
                    r_Archimedes_Spiral = exp_coefficients(1)*exp(-exp_coefficients(2)*t) + ...
                        exp_coefficients(3)*exp(-exp_coefficients(4)*t);
                else % i.e. for full optimization
                    r_Archimedes_Spiral = logspace(log10(r0),log10(rf),num_iter+1);% r_Archimedes_Spiral(end) = rf;
%                     if omega_init ~= 0
%                         r_Archimedes_Spiral = linspace(r0,rf,num_iter+1);
%                     else
%                         r_Archimedes_Spiral = [r0,r0-1e-10,rf+1e-10,rf];
%                         num_iter = 3;
%                     end
                end
                
                % Initialize variables
                Alin = zeros(3,1); Acor = zeros(3,1); Aang = zeros(3,1);
                Acen = zeros(3,1); Atot = zeros(3,1);
                A_TOT = zeros(1,1); A_LIN = zeros(1,1); A_COR = zeros(1,1);
                A_ANG = zeros(1,1); A_CEN = zeros(1,1);
                A_LIN_RATIO = zeros(1,1); A_COR_RATIO = zeros(1,1);
                A_ANG_RATIO = zeros(1,1); A_CEN_RATIO = zeros(1,1);
                DV_angular = zeros(1,1); DV_centripetal = zeros(1,1);
                DV_linear = zeros(1,1); DV_coriolis = zeros(1,1);
                DV_total = zeros(1,1);
                
                omega_init_norm = norm(omega_TAR_INT);
                
                for i = 1:num_iter
                    dt = dt_list(i);
                    q_UDP_INT(:,i+1) = mult_quat(q_UDP_TAR, q_TAR_INT(:,i)); %quaternion from inertial to UDP frame
                    q_SVC_INT(:,i+1) = conj_quat(q_UDP_INT(:,i+1));   %quaternion from inertial to chaser frame
                    
                    M_TAR(:,i+1) = [0;0;0];
                    tau_TAR(:,i+1) = [0;0;0];
                    
                    omega_dot_TAR_INT(:,i+1) = fast_dot(inv(I_TAR), (fast_cross(-omega_TAR_INT(:,i),I_TAR*omega_TAR_INT(:,i) + M_TAR(:,i+1)) + tau_TAR(:,i+1)));
                   
                    q_dot_TAR_INT(:,i+1) = 0.5 * mult_quat([omega_TAR_INT(:,i);0],q_TAR_INT(:,i));
                    
                    omega_TAR_INT(:,i+1) = omega_TAR_INT(:,i) + omega_dot_TAR_INT(:,i+1)*dt;
                    
                    quat_term = mult_quat(dt*q_dot_TAR_INT(:,i+1),conj_quat(q_TAR_INT(:,i)));
                    quat_term_vec_norm = norm(quat_term(1:3));
                    quat_exponential_term(:,i+1) = exp(quat_term(4)) * [quat_term(1:3)*sin(quat_term_vec_norm)/quat_term_vec_norm; cos(quat_term_vec_norm)];
                    q_TAR_INT(:,i+1) = mult_quat(quat_exponential_term(:,i+1),q_TAR_INT(:,i));
                    
                    v_SVC_INT(:,i+1) = fast_cross(omega_TAR_INT(:,i+1),r_SVC_INT(:,i));
                    
%                 if r_Archimedes_Spiral(i+1) < rf
%                     r_Archimedes_Spiral(i+1) = rf;
%                 end
                    r_SVC_INT(:,i+1) = trans_vec(q_TAR_INT(:,i+1), r_SVC_INT(:,1));% synchronous requirement!
                      

                    if work_with_exponentials == 0
                        r_SVC_INT(:,i+1) = r_SVC_INT(:,i+1) * r_Archimedes_Spiral(i+1)/norm(r_SVC_INT(:,i+1));
                        r_SVC_TAR(:,i+1) = [r_Archimedes_Spiral(:,i+1); 0; 0];
                        
                        % Matched Acceleration (2nd order differential
                        % equation)
%                         w1 = omega_TAR_INT(1,i+1); w2 = omega_TAR_INT(2,i+1); w3 = omega_TAR_INT(3,i+1);
%                         w1dot = omega_dot_TAR_INT(1,i+1); w2dot = omega_dot_TAR_INT(2,i+1); w3dot = omega_dot_TAR_INT(3,i+1);
%                         assignin('base','omega_TAR_INT',omega_TAR_INT)
%                         assignin('base','omega_dot_TAR_INT',omega_dot_TAR_INT)
%                         B = 2*[0 -w3 w2;...
%                             w3 0 -w1;...
%                             -w2 w1 0];
%                         C = [-w2^2-w3^2 -w3dot+w1*w2 w2dot+w1*w3;...
%                             w3dot+w1*w2 -w1^2-w3^2 -w1dot+w2*w3;...
%                             -w2dot+w1*w3 w1dot+w2*w3 -w2^2-w1^2];
%                         
%                         A = [zeros(3), eye(3);...
%                             C B];
%                         rdot_SVC_INT(:,i) = r_SVC_INT(:,i+1)-r_SVC_INT(:,i);
%                         x = [r_SVC_INT(:,i);rdot_SVC_INT(:,i)];
%                         dxdt(:,i) = A*x;
%                         assignin('base','dxdt',dxdt);
%                         dr(:,i) = dxdt(1:3,i);
%                         drdot(:,i) = dxdt(4:6,i);

                    else
                        r_SVC_INT(:,i+1) = r_SVC_INT(:,i+1) * r_Archimedes_Spiral/norm(r_SVC_INT(:,i+1));
                        r_SVC_TAR(:,i+1) = [r_Archimedes_Spiral; 0; 0];
                    end
    
                    rdot_SVC_TAR(:,i+1) = (r_SVC_TAR(:,i+1) - r_SVC_TAR(:,i)) /dt; 
                    % Calculate accelerations
                   if i < 2
%                         Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i)) /dt;
%                         A_LIN(i+1) = norm(Alin(:,i+1));
                        Acor(:,i+1) = 2* fast_cross(omega_TAR_INT(:,i+1),rdot_SVC_TAR(:,i+1));
                        A_COR(i+1) = norm(Acor(:,i+1));
                        Aang(:,i+1) = fast_cross(omega_dot_TAR_INT(:,i+1),r_SVC_TAR(:,i+1));
                        A_ANG(i+1) = norm(Aang(:,i+1));
                        Acen(:,i+1) = fast_cross(omega_TAR_INT(:,i+1), fast_cross(omega_TAR_INT(:,i+1),r_SVC_TAR(:,i+1)));
                        A_CEN(i+1) = norm(Acen(:,i+1));
                    else
%                         Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - 2*rdot_SVC_TAR(:,i) + rdot_SVC_TAR(:,i-1)) /dt^2;
%                         A_LIN(i+1) = norm(Alin(:,i+1));
                        Acor(:,i+1) = 2*fast_cross(omega_TAR_INT(:,i+1),rdot_SVC_TAR(:,i+1));
                        A_COR(i+1) = norm(Acor(:,i+1));
                        Aang(:,i+1) = fast_cross(omega_dot_TAR_INT(:,i+1),r_SVC_TAR(:,i+1));
                        A_ANG(i+1) = norm(Aang(:,i+1));
                        Acen(:,i+1) = fast_cross(omega_TAR_INT(:,i+1), fast_cross(omega_TAR_INT(:,i+1),r_SVC_TAR(:,i+1)));
                        A_CEN(i+1) = norm(Acen(:,i+1));
                    end
%                     if i == 1
%                         Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i)) /dt;
%                     elseif i == 2
%                         Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i-1)) /(2*dt);
%                     else
%                         Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i-2)) /(3*dt);
%                     end
                    Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i)) /dt;
                    A_LIN(i+1) = norm(Alin(:,i+1));
                    
                    Atot(:,i+1) = Alin(:,i+1) + Acor(:,i+1) + Aang(:,i+1) + Acen(:,i+1);
                    A_TOT(i+1) = sqrt(Atot(1,i+1)^2 + Atot(2,i+1)^2 + Atot(3,i+1)^2);
                    
                    A_LIN_RATIO(i+1) = A_LIN(i+1)/A_TOT(i+1);
                    A_COR_RATIO(i+1) = A_COR(i+1)/A_TOT(i+1);
                    A_ANG_RATIO(i+1) = A_ANG(i+1)/A_TOT(i+1);
                    A_CEN_RATIO(i+1) = A_CEN(i+1)/A_TOT(i+1);
                    
                    DV_angular(i+1)     = DV_angular(i) + A_ANG(:,i+1)*dt;
                    DV_centripetal(i+1) = DV_centripetal(i) + A_CEN(i+1)*dt;
                    DV_linear(i+1)      = DV_linear(i) + A_LIN(:,i+1)*dt;
                    DV_coriolis(i+1)    = DV_coriolis(i) + A_COR(i+1)*dt;
                    DV_total(i+1) = DV_total(i) + A_TOT(i+1)*dt;
                    
                    
                    %                 if r_Archimedes_Spiral(i+1) <= rf
                    %                     iterations = i;
                    %                     break;
                    %                 end
                    t(i+1) = t(i) + dt;
                end % end loop of simulating time until end of spiral
%                 DV_total = smooth(t,DV_total,0.1,'rloess');
                assignin('base','funcT',t)
                assignin('base','funcR',r_Archimedes_Spiral);
                
                if ~exist('iterations')
                    iterations = num_iter;
                end
                
                %             Alin(:,iterations+1:end) = [];
                %             Acor(:,iterations+1:end) = [];
                %             Aang(:,iterations+1:end) = [];
                %             Acen(:,iterations+1:end) = [];
                %             Atot(:,iterations+1:end) = [];
                %             A_TOT(:,iterations+1:end) = [];
                %             A_LIN(:,iterations+1:end) = [];
                %             A_COR(:,iterations+1:end) = [];
                %             A_ANG(:,iterations+1:end) = [];
                %             A_CEN(:,iterations+1:end) = [];
                %             DV_angular(:,iterations+1:end)      = [];
                %             DV_centripetal(:,iterations+1:end)  = [];
                %             DV_linear(:,iterations+1:end)       = [];
                %             DV_coriolis(:,iterations+1:end)     = [];
                %             DV_total(:,iterations+1:end)        = [];
                
                t_length = length(t);
                tf = t(end);
                
                DV_TOT = DV_total(end); DV_ANG = DV_angular(end);
                DV_CEN = DV_centripetal(end); DV_LIN = DV_linear(end);
                DV_COR = DV_coriolis(end);
                
                if length(A_TOT)>3
                    maxAtot = max(A_TOT);
                    AmaxIdx = find(A_TOT == maxAtot);
                    maxAlin = A_LIN(AmaxIdx); maxAcor = A_COR(AmaxIdx);
                    maxAcen = A_CEN(AmaxIdx); maxAang = A_ANG(AmaxIdx);
                else
                    maxAtot = A_TOT(end); maxAlin = A_LIN(end);
                    maxAcor = A_COR(end); maxAcen = A_CEN(end);
                    maxAang = A_ANG(end);
                end
            end %inertia loop
        end % rf_iter loop
    end % r0_iter loop
end % omega_iter loop

DeltaV = DV_TOT;

end