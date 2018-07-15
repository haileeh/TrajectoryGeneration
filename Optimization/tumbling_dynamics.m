q_UDP_INT(:,i+1) = mult_quat(q_UDP_TAR, q_TAR_INT(:,i)); %quaternion from inertial to UDP frame
q_SVC_INT(:,i+1) = conj_quat(q_UDP_INT(:,i+1));   %quaternion from inertial to chaser frame

% Create sine rotation multiply by regular omega
%omega_TAR_INT(:,i) =omega_TAR_INT(:,1).*(1+0.5*sin(5*t(i))); %omega_TAR_INT(:,i)+deg2rad(omega)*sin(t(i)); %since it's constant anyway - wont work for out of plane rotations

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

% Does this enforce synchrony?
r_SVC_INT(:,i+1) = trans_vec(q_TAR_INT(:,i+1), r_SVC_INT(:,1));   %does this make sense?
r_SVC_INT(:,i+1) = r_SVC_INT(:,i+1) * rad_dist(i+1)/norm(r_SVC_INT(:,i+1)); % This mult logspace distance by desired direction

r_SVC_TAR(:,i+1) = [rad_dist(:,i+1); 0; 0]; %maybe wrong
%r_SVC_TAR(:,i+1) = trans_vec(q_TAR_INT(:,i+1), r_SVC_INT(:,i+1));

rdot_SVC_TAR(:,i+1) = (r_SVC_TAR(:,i+1) - r_SVC_TAR(:,i)) /dt;