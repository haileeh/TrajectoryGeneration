function DeltaV = DeltaV_dt_r_omega0(dt_list,rad_dist,omega0, I_TAR_vec, plot_on_switch)
if nargin()<5
    plot_on = 1;
else
    plot_on = plot_on_switch;
end

rf = rad_dist(end); % m
r0 = rad_dist(1); % m

%Clear the data vectors
r_SVC_INT = [r0; 0; 0];  %m
r_SVC_TAR = r_SVC_INT;  % ASSUMPTION - keep it for now
r0 = norm(r_SVC_INT);
I_TAR = diag(I_TAR_vec);
omega_TAR_INT = deg2rad(omega0); %rate of target in inertial frame
omega = omega0;
v_SVC_INT = fast_cross(omega_TAR_INT,r_SVC_INT); %linear velocity of servicer relative to static target in inertial frame

q_TAR_INT = [0; 0; 0; 1];   %quaternion from inertial frame to target frame
q_UDP_TAR = [0; 0; 0; 1];   %quaternion from target frame to inertial frame

% Initialize variables
initialize_zeros
%rdot_SVC_TAR = [0.2 0 0]';%[0.01 0 0]';
num_iter = length(dt_list);
r = r0;
for i = 1:num_iter
    if r(i) < rf
        break;
    end
    dt = dt_list(i);
    tumbling_dynamics
    
    % Calculate accelerations and DV
    acc_DV_calculations
    
    % Normalize r to compare against rf
    r(i+1) = norm(r_SVC_TAR(:,i+1));
    
    assignin('base','Aang',Aang);
    assignin('base','Acor',Acor);
    assignin('base','Alin',Alin);
    assignin('base','Acen',Acen);
    assignin('base','Atot',Atot);
    
    t(i+1) = t(i) + dt;
end %end loop of simulating time until end of spiral

assignin('base','funcT',t)
assignin('base','funcR',rad_dist);

if ~exist('iterations')
    iterations = num_iter;
end

t_length = length(t);
tf = t(end);

DV_TOT = DV_total(end);
DV_ANG = DV_angular(end);
DV_CEN = DV_centripetal(end);
DV_LIN = DV_linear(end);
DV_COR = DV_coriolis(end);


%%
if plot_on == 1
    
    figure();
    plot(t(1:end),A_TOT(1:end),'k'); hold on;
    plot(t(1:end),A_LIN(1:end),'b');
    plot(t(1:end),A_COR(1:end),'r');
    plot(t(1:end),A_ANG(1:end),'m--');
    [hAx,hLine1,hLine2] = plotyy(t,A_CEN,t,norm(r_SVC_TAR));
    hLine1.LineStyle = ':';
    legend('Total Acceleration Magnitude','Linear','Coriolis','Angular','Centripetal');
    xlabel('t [s]')
    ylabel('Acceleration Component [m/s^2]')
    ylabel(hAx(2),'Radius [m]') % right y-axis
    grid on;
    
    figure();
    subplot(4,1,1)
    plot(t(1:end),A_TOT(1:end),'k'); hold on;
    plot(t(1:end),A_LIN(1:end)+A_CEN(1:end),'b');
    plot(t(1:end),A_COR(1:end)+A_ANG(1:end),'r');
    legend('Total','Linear+Centripetal','Coriolis+Angular');
    xlabel('t [s]')
    ylabel('Acceleration Norm')
    grid on;
    ylim([0,1.5]);
    subplot(4,1,2)
    plot(t(1:end),Alin(1,1:end)+Acen(1,1:end),'k'); hold on;
    plot(t(1:end),Acor(1,1:end)+Aang(1,1:end),'-.k');
    %plot(t(1:end),A_TOT(1:end),'g');
    legend('Lin1+Cen1','Cor1+Ang1');
    xlabel('t [s]')
    ylabel('Acceleration Component [m/s^2]')
    grid on;
    subplot(4,1,3)
    plot(t(1:end),Alin(2,1:end)+Acen(2,1:end),'b'); hold on;
    plot(t(1:end),Acor(2,1:end)+Aang(2,1:end),'-.b');
    legend('Lin2+Cen2','Cor2+Ang2');
    xlabel('t [s]')
    ylabel('Acceleration Component [m/s^2]')
    grid on;
    subplot(4,1,4)
    plot(t(1:end),Alin(3,1:end)+Acen(3,1:end),'r'); hold on;
    plot(t(1:end),Acor(3,1:end)+Aang(3,1:end),'-.r');
    legend('Lin3+Cen3','Cor3+Ang3');
    xlabel('t [s]')
    ylabel('Acceleration Component [m/s^2]')
    grid on;
    
    
    figure();
    angular_component = A_ANG(1:end)+A_COR(1:end);
    assignin('base','t',t)
    smoothed_radial_component = smooth(A_LIN(1:end) + A_CEN(1:end))';
    smoothed_radial_component(smoothed_radial_component>1) = 1+(1-angular_component(smoothed_radial_component>1));
    plot(t(1:end),smoothed_radial_component./(smoothed_radial_component+angular_component),'b'); hold on;
    [hAx,hLine1,hLine2] = plotyy(t,angular_component./(smoothed_radial_component+angular_component),t,r_SVC_TAR);
    hLine1.LineStyle = ':';
    legend('Radial','Tangential');
    xlabel('t [s]')
    ylabel({['Radial vs. Tangential Fraction']; ['of Total Acceleration']})
    ylabel(hAx(2),'Radius [m]') % right y-axis
    grid on;
    ylim([0 1.1])
    set(hAx(1),'XLim',[0 tf])
    set(hAx(2),'XLim',[0 tf])
      
    A_total_for_frac_smoothed = smooth(A_LIN+A_COR+A_ANG+A_CEN)';
    A_LIN_frac = A_LIN./A_total_for_frac_smoothed;
    A_COR_frac = A_COR./A_total_for_frac_smoothed;
    A_ANG_frac = A_ANG./A_total_for_frac_smoothed;
    A_CEN_frac = A_CEN./A_total_for_frac_smoothed;
    
    figure();
    plot(t(1:end),A_LIN_frac(1:end),'b'); hold on;
    plot(t(1:end),A_COR_frac(1:end),'r');
    plot(t(1:end),A_ANG_frac(1:end),'m--');
    plot(t,A_CEN_frac,':');
    legend('Linear','Coriolis','Angular','Centripetal');
    xlabel('t [s]')
    ylabel('Smoothed Acceleration Component Fraction of Total')
    grid on;
    ylim([0 1.1])
    xlim([0 tf])
    
    figure()
    plot(t(1:end),(A_LIN_frac+A_CEN_frac),'b'); hold on;
    [hAx,hLine1,hLine2] = plotyy(t,A_COR_frac+A_ANG_frac,t,norm(r_SVC_TAR));
    hLine1.LineStyle = ':';
    legend('Radial','Tangential');
    xlabel('t [s]')
    ylabel({['Radial vs. Tangential Fraction']; ['of Total Acceleration']})
    ylabel(hAx(2),'Radius [m]') % right y-axis
    grid on;
    ylim([0 1.1])
    set(hAx(1),'XLim',[0 tf])
    set(hAx(2),'XLim',[0 tf])
    
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
    xlabel('t [s]')
    ylabel('Acceleration Component Fraction of Total')
    grid on;
    ylim([0 1.1])
    xlim([0 tf])
   
    for i = 1:iterations
        [yaw(i), pitch(i), roll(i)] = euler_from_quat(q_TAR_INT(:,i));
        dcm = dcm_from_quat(q_TAR_INT(:,i));
        DCM(i,:,:) = dcm;
        pointer_TAR_INT(:,i) = dcm*[1;0;0];
    end
    
    figure();
    plot3(r_SVC_INT(1,2:iterations),r_SVC_INT(2,2:iterations),r_SVC_INT(3,2:iterations)); hold on;
    plot3(pointer_TAR_INT(1,1:iterations),pointer_TAR_INT(2,1:iterations),pointer_TAR_INT(3,1:iterations));
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
    
end

assignin('base','q_TAR_INT',q_TAR_INT)
assignin('base','omega_TAR_INT',omega_TAR_INT)
assignin('base','omega_dot_TAR_INT',omega_dot_TAR_INT)
DeltaV = DV_TOT;

end
