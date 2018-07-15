function DeltaV = spiralDeltaV_dt_MultipleModel_exponentials(exponential_parameters,omega,R0,RF,inertia_ratios,tf,num_iter)

% Initialize Output
DeltaV = 0;

OMEGA = size(omega,2); %deg/sec
INERTIA = size(inertia_ratios,2);

dt_list = tf/num_iter * ones(num_iter,1);
t_list = cumsum(dt_list)';
if length(exponential_parameters) == 8
    r_Archimedes_Spiral = exponential_parameters(1)*exp(-exponential_parameters(2)*t_list) + ...
        exponential_parameters(3)*exp(-exponential_parameters(4)*t_list) + ...
        exponential_parameters(5)*exp(-exponential_parameters(6)*t_list) + ...
        exponential_parameters(7)*exp(-exponential_parameters(8)*t_list);
end
for omega_iter = 1:(OMEGA)
    for r0_iter = 1:length(R0)
        for rf_iter = 1:length(RF)
            for INERTIA_idx = 1:(INERTIA)
                
                %Set the trajectory parameters
                rf = RF(rf_iter);
                r0 = R0(r0_iter);
                omega_init = norm(omega(:,omega_iter));
                
                %Clear the data vectors
                r_SVC_INT = [r0; 0; 0];  %m
                r_SVC_TAR = r_SVC_INT;
                r0 = norm(r_SVC_INT);
                
                I_TAR = diag(inertia_ratios(:,INERTIA_idx));
                omega_TAR_INT = deg2rad(omega); %rate of target in inertial frame
                
                
                v_SVC_INT = fast_cross(omega_TAR_INT,r_SVC_INT); %linear velocity of servicer relative to static target in inertial frame
                q_TAR_INT = [0; 0; 0; 1];   %quaternion from inertial frame to target frame
                q_UDP_TAR = [0; 0; 0; 1];   %quaternion from target frame to inertial frame
%                 q_UDP_TAR = [0.0019;   -0.0436;    0.0436;    0.9981];
%                 q_UDP_TAR = [0 0.7071 0.7071 0]';
                %                                 num_iter = 1000000;
%                 num_iter = length(dt_list);
                q_UDP_INT = zeros(4,1);
                q_SVC_INT = zeros(4,1);
                M_TAR = zeros(3,1);
                tau_TAR = zeros(3,1);
                omega_dot_TAR_INT = zeros(3,1);
                q_dot_TAR_INT = zeros(4,1);
                r_TEMP = zeros(3,1);
%                 r_Archimedes_Spiral = linspace(r0,rf,num_iter+1);
%                 r_Archimedes_Spiral = logspace(log10(r0),log10(rf),num_iter+1);% r_Archimedes_Spiral(end) = rf;
                t = 0;
                Alin = zeros(3,1);
                Acor = zeros(3,1);
                Aang = zeros(3,1);
                Acen = zeros(3,1);
                Atot = zeros(3,1);
                A_TOT = zeros(1,1);
                A_LIN = zeros(1,1);
                A_COR = zeros(1,1);
                A_ANG = zeros(1,1);
                A_CEN = zeros(1,1);
                A_LIN_RATIO = zeros(1,1);
                A_COR_RATIO = zeros(1,1);
                A_ANG_RATIO = zeros(1,1);
                A_CEN_RATIO = zeros(1,1);
                DV_angular = zeros(1,1);
                DV_centripetal = zeros(1,1);
                DV_linear = zeros(1,1);
                DV_coriolis = zeros(1,1);
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
                    r_SVC_INT(:,i+1) = trans_vec(q_TAR_INT(:,i+1), r_SVC_INT(:,1));
                    r_SVC_INT(:,i+1) = r_SVC_INT(:,i+1) * r_Archimedes_Spiral(i)/norm(r_SVC_INT(:,i+1));
                    r_SVC_TAR(:,i+1) = [r_Archimedes_Spiral(:,i); 0; 0];
                    
                    rdot_SVC_TAR(:,i+1) = (r_SVC_TAR(:,i+1) - r_SVC_TAR(:,i)) /dt;
                    
                    
                    % Acceleration Terms: Inertial Frame
                    Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i)) /dt;
                    A_LIN(i+1) = norm(Alin(:,i+1));
                    Acor(:,i+1) = 2* fast_cross(omega_TAR_INT(:,i+1),rdot_SVC_TAR(:,i+1));
                    A_COR(i+1) = norm(Acor(:,i+1));
                    Aang(:,i+1) = fast_cross(omega_dot_TAR_INT(:,i+1),r_SVC_TAR(:,i+1));
                    A_ANG(i+1) = norm(Aang(:,i+1));
                    Acen(:,i+1) = fast_cross(omega_TAR_INT(:,i+1), fast_cross(omega_TAR_INT(:,i+1),r_SVC_TAR(:,i+1)));
                    A_CEN(i+1) = norm(Acen(:,i+1));
                    
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
                    
                    
                    t(i+1) = t(i) + dt;
                end %end loop of simulating time until end of spiral
                
                assignin('base','funcT',t)
                assignin('base','funcR',r_Archimedes_Spiral);
                
                if ~exist('iterations')
                    iterations = num_iter;
                end
                
                
                t_length = length(t);
                tf = t(end);
                
%                 DV_total = smooth(t,DV_total,0.1,'rloess');
                DV_TOT = DV_total(end);
                DV_ANG = DV_angular(end);
                DV_CEN = DV_centripetal(end);
                DV_LIN = DV_linear(end);
                DV_COR = DV_coriolis(end);
                
                if length(A_TOT)>3
                    maxAtot = max(A_TOT);
                    AmaxIdx = find(A_TOT == maxAtot);
                    maxAlin = A_LIN(AmaxIdx);
                    maxAcor = A_COR(AmaxIdx);
                    maxAcen = A_CEN(AmaxIdx);
                    maxAang = A_ANG(AmaxIdx);
                else
                    maxAtot = A_TOT(end);
                    maxAlin = A_LIN(end);
                    maxAcor = A_COR(end);
                    maxAcen = A_CEN(end);
                    maxAang = A_ANG(end);
                end
                
                DeltaV = DeltaV + DV_TOT;
                
            end %inertia loop
        end % rf_iter loop
    end % r0_iter loop
end % omega_iter loop

DeltaV = DeltaV / ((OMEGA) * length(R0) * length(RF) * (INERTIA));

end