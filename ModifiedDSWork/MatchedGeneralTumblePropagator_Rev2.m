clear pointer_TAR_INT
clear dr
clear drdot
clc; close all;

%% Setup
num_iter = 2000; %100000;
dt = 0.1; % sec

% quaternion setup
% q_TAR_INT is only initial - should update
q_TAR_INT = [0; 0; 0; 1];   %quaternion from inertial frame to target frame
% Above assumes that TAR = INT initially (which is unlikely)
% q_UDP_TAR is fixed!
q_UDP_TAR = [0; 0; 0; 1];   %quaternion from UDP frame to Target frame
%q_UDP_TAR = [0.0019;   -0.0436;    0.0436;    0.9981];
%q_UDP_TAR = [0 0.7071 0.7071 0]';

r0 = 10;
r = r0;
rf = 1;
r_SVC_INT = [r0; 0; 0];  %[r0/sqrt(3);r0/sqrt(3);r0/sqrt(3)]; %[r0; 0; 0]; 

%% Different rotations setup
inertia_choice = 'fully_symmetric';
rot_choice = 'w3_rot';
% Principle Axis Rotations
switch inertia_choice
    case 'fully_symmetric'  % always correlates to single axis rotation?
       I_TAR = diag([1 1 1]'); 
    case 'axis_sym_1'
        I_TAR = diag([1.8534 1 1]');
        %I_TAR = diag([1.95 1 1]');
        %I_TAR = diag([1.1990 1 1]');
    case 'axis_sym_3'
        I_TAR = diag([1.8534 1.8534 1]'); 
        %I_TAR = diag([2 2 1]'); 
    case 'tri_sym'
        I_TAR = diag([1.1990 1.1642 1]');
end

switch rot_choice
    case 'w1_rot'
        omega_TAR_INT=deg2rad( 5/norm([1 0 0])*[1 0 0]');
        r_SVC_INT = [0; r0; 0];
    case 'w2_rot'
        omega_TAR_INT=deg2rad( 5/norm([0 1 0])*[0 1 0]');
    case 'w3_rot'
        omega_TAR_INT=deg2rad( 5/norm([0 0 1])*[0 0 1]');
    case 'multi_rot1'
        omega_TAR_INT = deg2rad( 5/norm([0.5 10 0.5])*[0.5 10 0.5]'); %rad/s
    case 'multi_rot2'
        omega_TAR_INT = deg2rad( 5/norm([0 10 0.5])*[0 10 0.5]'); %rad/s
    case 'multi_rot3'
        omega_TAR_INT = deg2rad( 5/norm([10 0 0.5])*[10 0 0.5]'); %rad/s
        r_SVC_INT = [0; r0; 0];
end

% check energy case
J = I_TAR/I_TAR(3,3);
J1 = J(1,1); J2 = J(2,2); J3 = J(3,3);
w1 = omega_TAR_INT(1); w2 = omega_TAR_INT(2); w3 = omega_TAR_INT(3);
H = J1^2*w1^2 + J2^2*w2^2 + J3^2*w3^2;
E = J1*w1^2 + J2*w2^2 + J3*w3^2;
energy_case = H - E*J2; % 0 = medium; + = low; - = high

r_SVC_TAR = trans_vec(q_TAR_INT,r_SVC_INT);
%r_SVC_TAR = r_SVC_INT; % assumes TAR and INT are the same frame?

initialize_zeros
% backward prop = use negative step size and define rdot (final 

for i = 1:num_iter
    if r(i) < rf
        break;
    end
    
    q_UDP_INT(:,i+1) = mult_quat(q_UDP_TAR, q_TAR_INT(:,i)); %quaternion from inertial to UDP frame
    % supposed to be rotating it about z-axis to obtain the quaternion the
    % servicer must obtain for both docking ports to align!
    q_SVC_INT(:,i+1) = conj_quat(q_UDP_INT(:,i+1));   %quaternion from inertial to chaser frame
    
    M_TAR(:,i+1) = [0;0;0];
    tau_TAR(:,i+1) = [0;0;0];
    
    % won't change unless M or tau acts on it
    omega_dot_TAR_INT(:,i+1) = fast_dot(inv(I_TAR), (fast_cross(-omega_TAR_INT(:,i),I_TAR*omega_TAR_INT(:,i) + M_TAR(:,i+1)) + tau_TAR(:,i+1)));    
    omega_TAR_INT(:,i+1) = omega_TAR_INT(:,i) + omega_dot_TAR_INT(:,i+1)*dt;
    
    %  All to update q_TAR_INT
    q_dot_TAR_INT(:,i+1) = 0.5*mult_quat([omega_TAR_INT(:,i);0],q_TAR_INT(:,i)); % quaternion derivative
    quat_term = mult_quat(dt*q_dot_TAR_INT(:,i+1),conj_quat(q_TAR_INT(:,i)));
    % 3 lines below are definitely correct
    quat_term_vec_norm = norm(quat_term(1:3));
    quat_exponential_term(:,i+1) = exp(quat_term(4)) * [quat_term(1:3)*sin(quat_term_vec_norm)/quat_term_vec_norm; cos(quat_term_vec_norm)];
    q_TAR_INT(:,i+1) = mult_quat(quat_exponential_term(:,i+1),q_TAR_INT(:,i));
    
    % This should happen in Target's frame
    w1(i) = omega_TAR_INT(1,i); w2(i) = omega_TAR_INT(2,i); w3(i) = omega_TAR_INT(3,i);
    w1dot(i) = omega_dot_TAR_INT(1,i); w2dot(i) = omega_dot_TAR_INT(2,i); w3dot(i) = omega_dot_TAR_INT(3,i);
    
    w_skew = [0 -w3(i) w2(i) ; w3(i) 0 -w1(i) ; -w2(i) w1(i) 0];
    wdot_skew = [0 -w3dot(i) w2dot(i) ; w3dot(i) 0 -w1dot(i) ; -w2dot(i) w1dot(i) 0];
    B = 2*w_skew;
    C = (wdot_skew - w_skew'*w_skew)/dt; %added transpose, fixed signs
    A = [zeros(3), eye(3);...
        C B];
    x = [r_SVC_TAR(:,i);rdot_SVC_TAR(:,i)];
    dxdt(:,i) = A*x;
    dr(:,i) = dxdt(1:3,i);
    drdot(:,i) = dxdt(4:6,i);
    %

    r_SVC_TAR(:,i+1) = r_SVC_TAR(:,i) + drdot(:,i)*(1*dt);
    
    %r_SVC_TAR(1,i+1) = r_SVC_TAR(1,i) + drdot(1,i)*(1*dt);  % Sternberg
    
    rdot_SVC_TAR(:,i+1) = (r_SVC_TAR(:,i+1) - r_SVC_TAR(:,i)) / (1*dt);
    
    %v_SVC_INT(:,i+1) = fast_cross(omega_TAR_INT(:,i+1),r_SVC_INT(:,i));
    
    % what does below even do?
    %r_SVC_INT(:,i+1) = trans_vec(q_TAR_INT(:,i+1), r_SVC_INT(:,1));% synchronous requirement?
    r_SVC_INT(:,i+1) = trans_vec(conj_quat(q_TAR_INT(:,i+1)), r_SVC_TAR(:,i+1));
    % Assuming that TAR and UDP are coaligned. This rotates vector into the TAR frame
    % without this req, need to actively command r_SVC_INT??
    
    % Normalize r to compare against rf
    r(i+1) = norm(r_SVC_TAR(:,i+1));
    
    acc_DV_calculations
    
    t(i+1) = t(i) + dt;
end %end loop of simulating time until end of spiral

t_length = length(t);
tf = t(end);

DV_TOT = DV_total(end)
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

figure(); plot(t,r); hold on; grid on;
xlabel('Time [s]')
ylabel('Radius [m]')

A_total_for_frac = A_LIN(1:end)+A_COR(1:end)+A_ANG(1:end)+A_CEN(1:end);
A_LIN_frac = A_LIN(1:end)./A_total_for_frac;
A_COR_frac = A_COR(1:end)./A_total_for_frac;
A_ANG_frac = A_ANG(1:end)./A_total_for_frac;
A_CEN_frac = A_CEN./A_total_for_frac;

figure();
plot(t(1:end),A_LIN_frac(1:end),'b'); hold on;
plot(t(1:end),A_COR_frac(1:end),'r');
plot(t(1:end),A_ANG_frac(1:end),'m--');
plot(t,A_CEN_frac,':');
legend('Linear','Coriolis','Angular','Centripetal');
xlabel('Time [s]')
ylabel('Acceleration Component Fraction of Total')
grid on;
ylim([0 1.1])
%xlim([0 tf])

figure();
plot(t(1:end),A_LIN_frac(1:end)+A_CEN_frac,'b'); hold on;
plot(t(1:end),A_ANG_frac+A_COR_frac(1:end),'r');
legend('Radial','Tangential');
xlabel('t [s]')
ylabel('Acceleration Component Fraction of Total')
grid on;
ylim([0 1.1])
%xlim([0 tf])

for i = 1:t_length
    [yaw, pitch, roll] = euler_from_quat(q_TAR_INT(:,i));
    dcm = dcm_from_quat(q_TAR_INT(:,i));
    % Illustrates wobble? how TAR frame changes wrt INT frame?
    pointer_TAR_INT(:,i) = dcm*[1;0;0]; % selects first column?
end

figure
plot3(r_SVC_INT(1,2:end),r_SVC_INT(2,2:end),r_SVC_INT(3,2:end)); hold on;
plot3(pointer_TAR_INT(1,1:end),pointer_TAR_INT(2,1:end),pointer_TAR_INT(3,1:end));
[X,Y,Z] = sphere(50);
X = 0.99*X;
Y = 0.99*Y;
Z = 0.99*Z;
sh = surfl(X, Y, Z); % capture the handle, i.e. the unique identifier, of the surface
% set color to gray, make mostly transparent
set(sh,'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.75)
axis equal;
grid on;
title({['Chaser Approach to Target with'];['r_0=' , num2str(r0), '[m], and r_f=', num2str(rf), '[m]']})
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');

figure
plot(t,pointer_TAR_INT(1,1:end),t,pointer_TAR_INT(2,1:end),t,pointer_TAR_INT(3,1:end))
title('pointer TAR INT time history');
xlabel('Time [s]');
ylabel('Motion');
