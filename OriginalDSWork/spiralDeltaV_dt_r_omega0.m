function DeltaV = spiralDeltaV_dt_r_omega0(dt_list,r_list_in,omega0, I_TAR_vec, plot_on_switch)
if nargin()<5
    plot_on = 1;
else
    plot_on = plot_on_switch;
end
size(r_list_in)
r_list =[r_list_in, r_list_in(end)];

RF = 1; %m
R0 = 10;%m
OMEGA = norm(omega0); %deg/sec
INERTIA = size(I_TAR_vec,2);

for omega_iter = 1:length(OMEGA)
    for r0_iter = 1:length(R0)
        for rf_iter = 1:length(RF)
            
            %Set the trajectory parameters
            rf = RF(rf_iter);
            r0 = R0(r0_iter);
            omega_init = omega0;
            
            %Clear the data vectors
            r_SVC_INT = [r0; 0; 0];  %m
            r_SVC_TAR = r_SVC_INT;
            r0 = norm(r_SVC_INT);
            %             I_TAR = [1 0 0;...      %kgm^2
            %                 0 1 0;...
            %                 0 0 1];
            %                     I_TAR = [0.028856093583883000,-0.000271760426029421, 0.001061424114444470;...
            %                                               -0.000271760426029421, 0.055986602544449500, 0.000183704942897512;...
            %                                                0.001061424114444470, 0.000183704942897512, 0.058822298355749900];
            %                                 I_TAR = [3,-0.00271760426029421, 0.01061424114444470;...
            %                                     -0.00271760426029421, 1, 0.00183704942897512;...
            %                                     0.01061424114444470, 0.00183704942897512, 2];
            %                     I_TAR = diag([1.239 1.1905 1]); % Tim's original
            %             I_TAR = diag([3 2 1]);
            I_TAR = diag(I_TAR_vec);
            %                     omega_TAR_INT = deg2rad([omega_init;0;0]); %rate of target in inertial frame
            %                     omega_TAR_INT = deg2rad([0; omega_init;0]); %rate of target in inertial frame
            omega_TAR_INT = deg2rad(omega_init); %rate of target in inertial frame
            %                     omega_TAR_INT = [0.699809513728002 0.716087962216080 0.355816171813249]'; %Tim's original
            %                     omega_TAR_INT = [0.07; 0.072; 0.036]; %Rounded 0.1*Tim's original
            
            v_SVC_INT = fast_cross(omega_TAR_INT,r_SVC_INT); %linear velocity of servicer relative to static target in inertial frame
            q_TAR_INT = [0; 0; 0; 1];   %quatenrion from inertial frame to target frame
            %             q_UDP_TAR = [0; 0; 0; 1];   %quaternion from target frame to inertial frame
            %             q_UDP_TAR = [0.0019;   -0.0436;    0.0436;    0.9981];
            q_UDP_TAR = [0 0.7071 0.7071 0]';
            
            %                                 num_iter = 1000000;
            num_iter = length(dt_list);
            q_UDP_INT = zeros(4,1);
            q_SVC_INT = zeros(4,1);
            M_TAR = zeros(3,1);
            tau_TAR = zeros(3,1);
            omega_dot_TAR_INT = zeros(3,1);
            q_dot_TAR_INT = zeros(4,1);
            r_TEMP = zeros(3,1);
            r_Archimedes_Spiral = r_list;
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
                %                 q_TAR_INT(:,i+1)
                %                 r_SVC_INT(:,i+1)
                r_SVC_INT(:,i+1) = r_SVC_INT(:,i+1) * r_Archimedes_Spiral(i+1)/norm(r_SVC_INT(:,i+1));
                r_SVC_TAR(:,i+1) = [r_Archimedes_Spiral(:,i+1); 0; 0];
                
                rdot_SVC_TAR(:,i+1) = (r_SVC_TAR(:,i+1) - r_SVC_TAR(:,i)) /dt;
                
%                 Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i)) /dt;
%                 Acceleration Terms: Inertial Frame
                if i == 1
                    Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i)) /dt;
                elseif i == 2
                    Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i-1)) /(2*dt);
                else
                    Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i-2)) /(3*dt);
%                 else
%                     Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i-3)) /(4*dt);
                end
                %                 if norm(Alin(:,i+1)) > 10*norm(Alin(:,i))
                %                     Alin(:,i+1) = Alin(:,i);
                %                 end
                A_LIN(i+1) = norm(Alin(:,i+1));
                Acor(:,i+1) = 2* fast_cross(omega_TAR_INT(:,i+1),rdot_SVC_TAR(:,i+1));
                A_COR(i+1) = norm(Acor(:,i+1));
                Aang(:,i+1) = fast_cross(omega_dot_TAR_INT(:,i+1),r_SVC_TAR(:,i+1));
                A_ANG(i+1) = norm(Aang(:,i+1));
                Acen(:,i+1) = fast_cross(omega_TAR_INT(:,i+1), fast_cross(omega_TAR_INT(:,i+1),r_SVC_TAR(:,i+1)));
                A_CEN(i+1) = norm(Acen(:,i+1));
                
                Atot(:,i+1) = Alin(:,i+1) + Acor(:,i+1) + Aang(:,i+1) + Acen(:,i+1);
                A_TOT(i+1) = sqrt(Atot(1,i+1)^2 + Atot(2,i+1)^2 + Atot(3,i+1)^2);
                
                assignin('base','Aang',Aang);
                assignin('base','Acor',Acor);
                assignin('base','Alin',Alin);
                assignin('base','Acen',Acen);
                assignin('base','Atot',Atot);
                
                A_LIN_RATIO(i+1) = A_LIN(i+1)/A_TOT(i+1);
                A_COR_RATIO(i+1) = A_COR(i+1)/A_TOT(i+1);
                A_ANG_RATIO(i+1) = A_ANG(i+1)/A_TOT(i+1);
                A_CEN_RATIO(i+1) = A_CEN(i+1)/A_TOT(i+1);
                
                DV_angular(i+1)     = DV_angular(i) + A_ANG(:,i+1)*dt;
                DV_centripetal(i+1) = DV_centripetal(i) + A_CEN(i+1)*dt;
                DV_linear(i+1)      = DV_linear(i) + A_LIN(:,i+1)*dt;
                DV_coriolis(i+1)    = DV_coriolis(i) + A_COR(i+1)*dt;
                DV_total(i+1) = DV_total(i) + A_TOT(i+1)*dt;
%                 if norm(Alin(:,i+1)) > 10*norm(Alin(:,i)) && i>4
%                     Alin(:,i+1) = Alin(:,i);
%                 end
%                 A_LIN(i+1) = norm(Alin(:,i+1));
                
                %                 if r_Archimedes_Spiral(i+1) <= rf
                %                     iterations = i;
                %                     break;
                %                 end
                t(i+1) = t(i) + dt;
            end %end loop of simulating time until end of spiral
            
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
            
            assignin('base','dt_list',dt_list);
            assignin('base','rdot_SVC_TAR',rdot_SVC_TAR);
            assignin('base','r_SVC_TAR',r_SVC_TAR);
            assignin('base','r_Archimedes_Spiral',r_Archimedes_Spiral);
            assignin('base','r_SVC_INT',r_SVC_INT);
            assignin('base','v_SVC_INT',v_SVC_INT);
            
        end % rf_iter loop
    end % r0_iter loop
end % omega_iter loop
%%
if plot_on == 1    
    figure();
    plot(t(1:end),A_TOT(1:end),'k'); hold on;
    plot(t(1:end),A_LIN(1:end),'b');
    plot(t(1:end),A_COR(1:end),'r');
    plot(t(1:end),A_ANG(1:end),'m--');
    [hAx,hLine1,hLine2] = plotyy(t,A_CEN,t,r_SVC_TAR);
    hLine1.LineStyle = ':';
    legend('Total Acceleration Magnitude','Linear','Coriolis','Angular','Centripetal');
    xlabel('t [s]')
    ylabel('Acceleration Component [m/s^2]')
    ylabel(hAx(2),'Radius [m]') % right y-axis
    grid on;
    
    figure();
    subplot(2,1,1)
    plot(t(1:end),A_TOT(1:end),'k'); hold on;
    plot(t(1:end),A_LIN(1:end)+A_CEN(1:end),'b');
    plot(t(1:end),A_COR(1:end),'r');
    legend('Total','Linear+Centripetal','Coriolis');
    xlabel('t [s]')
    ylabel('Acceleration Component [m/s^2]')
    grid on;
    ylim([0,1.5]);
    subplot(2,1,2)
    plot(t(1:end),Alin(1,1:end),'k'); hold on;
    plot(t(1:end),Alin(2,1:end),'b');
    plot(t(1:end),Alin(3,1:end),'r');
    plot(t(1:end),Acen(1,1:end),'--k');
    plot(t(1:end),Acen(2,1:end),'--b');
    plot(t(1:end),Acen(3,1:end),'--r');
    plot(t(1:end),Acor(1,1:end),'-.k');
    plot(t(1:end),Acor(2,1:end),'-.b');
    plot(t(1:end),Acor(3,1:end),'-.r');
    plot(t(1:end),A_TOT(1:end),'g');
    legend('Lin1','Lin2','Lin3','Cen1','Cen2','Cen3','Cor1','Cor2','Cor3','Total Acceleration Magnitude');
    xlabel('t [s]')
    ylabel('Acceleration Component [m/s^2]')
    grid on;
    
    figure();
    assignin('base','A_LIN',A_LIN)
    assignin('base','A_CEN',A_CEN)
    assignin('base','A_ANG',A_ANG)
    assignin('base','A_COR',A_COR)
    assignin('base','A_TOT',A_TOT)
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
    
    % This plot doesn't seem to provide a "sum to 1" if you add all the
    % components together. Need to check.
%     figure();
%     plot(t(1:end),A_LIN_RATIO(1:end),'b'); hold on;
%     plot(t(1:end),A_COR_RATIO(1:end),'r');
%     plot(t(1:end),A_ANG_RATIO(1:end),'m--');
%     [hAx,hLine1,hLine2] = plotyy(t,A_CEN_RATIO,t,r_SVC_TAR);
%     hLine1.LineStyle = ':';
%     legend('Linear','Coriolis','Angular','Centripetal');
%     xlabel('t [s]')
%     ylabel('Acceleration Component Fraction of Total')
%     ylabel(hAx(2),'Radius [m]') % right y-axis
%     grid on;
%     ylim([0 1.5])
%     set(hAx(1),'XLim',[0 tf])
%     set(hAx(2),'XLim',[0 tf])

    
    A_total_for_frac_smoothed = smooth(A_LIN+A_COR+A_ANG+A_CEN)';
    A_LIN_frac = A_LIN./A_total_for_frac_smoothed;
    A_COR_frac = A_COR./A_total_for_frac_smoothed;
    A_ANG_frac = A_ANG./A_total_for_frac_smoothed;
    A_CEN_frac = A_CEN./A_total_for_frac_smoothed;
    assignin('base','A_LIN_frac',A_LIN_frac)
    assignin('base','A_CEN_frac',A_CEN_frac)
    assignin('base','A_ANG_frac',A_ANG_frac)
    assignin('base','A_COR_frac',A_COR_frac)
    assignin('base','A_total_for_frac_smoothed',A_total_for_frac_smoothed')
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
    
    figure()
    plot(t(1:end),1-(A_COR_frac+A_ANG_frac),'b'); hold on;
    [hAx,hLine1,hLine2] = plotyy(t,A_COR_frac+A_ANG_frac,t,r_SVC_TAR);
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
    assignin('base','A_LIN_frac',A_LIN_frac)
    assignin('base','A_CEN_frac',A_CEN_frac)
    assignin('base','A_ANG_frac',A_ANG_frac)
    assignin('base','A_COR_frac',A_COR_frac)
    assignin('base','A_total_for_frac_smoothed',A_total_for_frac')
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

    
    figure();
    plot(t(2:end),A_LIN_RATIO(2:end)./A_LIN_RATIO(2:end),'b'); hold on;
    plot(t(2:end),A_COR_RATIO(2:end)./A_LIN_RATIO(2:end),'r');
    plot(t(2:end),A_ANG_RATIO(2:end)./A_LIN_RATIO(2:end),'m--');
    [hAx,hLine1,hLine2] = plotyy(t,A_CEN_RATIO./A_LIN_RATIO,t,r_SVC_TAR);
    hLine1.LineStyle = ':';
    legend('Linear','Coriolis','Angular','Centripetal');
    xlabel('t [s]')
    ylabel({['Acceleration Normalized by']; ['Linear Component']})
    ylabel(hAx(2),'Radius [m]') % right y-axis
    grid on;
    
    
    for i = 1:iterations
        [yaw, pitch, roll] = euler_from_quat(q_TAR_INT(:,i));
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
    
    figure();
    [hAx,hLine1,hLine2] = plotyy(t,DV_total,t,r_SVC_TAR);
    hLine1.LineStyle = ':';
    legend('\DeltaV','Radial Distance');
    xlabel('t [s]')
    ylabel('\DeltaV [m/s]')
    ylabel(hAx(2),'Radius [m]') % right y-axis
    grid on;
    
    figure();
    plot(t(1:end),A_TOT(1:end)./max(A_TOT),'k'); hold on;
    [hAx,hLine1,hLine2] = plotyy(t,DV_total./DV_total(end),t,r_SVC_TAR);
    hLine1.LineStyle = ':';
    legend('Normalized A_{TOT}','Normalized \DeltaV','Radial Distance');
    xlabel('t [s]')
    ylabel('Normalized \DeltaV and Total Acceleration')
    ylabel(hAx(2),'Radius [m]') % right y-axis
    grid on;
    
    % figure();
    % d = designfilt('lowpassfir','PassbandFrequency',0.25, ...
    %          'StopbandFrequency',0.35,'PassbandRipple',0.5, ...
    %          'StopbandAttenuation',65,'DesignMethod','kaiserwin');
    % t_new = interp1(linspace(0,t(end),length(t)),t,linspace(0,t(end),10*length(t)));
    % A_TOT = interp1(t,A_TOT,t_new);
    % A_LIN = interp1(t,A_LIN,t_new);
    % A_COR = interp1(t,A_COR,t_new);
    % A_ANG = interp1(t,A_ANG,t_new);
    % plot(t_new,filtfilt(d,A_TOT),'k'); hold on;
    % plot(t_new,filtfilt(d,A_LIN),'b');
    % plot(t_new,filtfilt(d,A_COR),'r');
    % plot(t_new,filtfilt(d,A_ANG),'m--');
    % [hAx,hLine1,hLine2] = plotyy(t,A_CEN,t,r_SVC_TAR);
    % hLine1.LineStyle = ':';
    % legend('Total','Linear','Coriolis','Angular','Centripetal');
    % xlabel('t [s]')
    % ylabel('Acceleration Component [m/s^2]')
    % ylabel(hAx(2),'Radius [m') % right y-axis
    % grid on;
end

assignin('base','q_TAR_INT',q_TAR_INT)
assignin('base','omega_TAR_INT',omega_TAR_INT)
assignin('base','omega_dot_TAR_INT',omega_dot_TAR_INT)
DeltaV = DV_TOT;

end