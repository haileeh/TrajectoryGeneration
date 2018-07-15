% NOTE: This script must be run as CELLS ONLY and AFTER
% Feasibility_Space_Archimedes_Spirals_Multiple_rf_and_omega.m

%% Find Feasible Targets for a Chaser
return_val = 1;
idx = linspace(1,size(TOTAL_POINTS,2),size(TOTAL_POINTS,2));
while return_val == 1
    % Plot to allow user to select which properties to plot and constrain
    figure(888)
    plot([1 2 3 4], [0 0 0 0], '*')
    title('Select Chaser Properties to Constrain')
    xlim([0 5])
    set(gca,'YTickLabel',{' '})
    set(gca,'XTickLabel',{'', '\DeltaV_{tot}', 'a_{max}', 't_f', 'v_f', ''});
    [ix,~] = ginput(2);
    ix=round(ix,0);
    total_points_index = ix + 2;
    close Figure 888
    
    % Plot enabling the constraint values to be picked
    figure(889)
    scatter(TOTAL_POINTS(total_points_index(1),:),TOTAL_POINTS(total_points_index(2),:));
    title({'Constrain Chaser Satellite Properties:'; 'Click 1 is for X Axis, Click 2 is for Y Axis'})
    if total_points_index(1) == 3
        xlabel('\DeltaV Required [m/s]');
    elseif total_points_index(1) == 4
        xlabel('Peak Acceleration [m/s^2]');
    elseif total_points_index(1) == 5
        xlabel('Final Time [s]');
    elseif total_points_index(1) == 6
        xlabel('Final Velocity [m/s]');
    end
    
    if total_points_index(2) == 3
        ylabel('\DeltaV Required [m/s]');
    elseif total_points_index(2) == 4
        ylabel('Peak Acceleration [m/s^2]');
    elseif total_points_index(2) == 5
        ylabel('Final Time [s]');
    elseif total_points_index(2) == 6
        ylabel('Final Velocity [m/s]');
    end
    grid on;
    [x,y] = ginput(2);
    close Figure 889
    
    [~,idx_this_pass] = find(TOTAL_POINTS(total_points_index(1),:) <= x(1) & TOTAL_POINTS(total_points_index(2),:) <= y(2));
    Lia = ismember(idx, idx_this_pass);
    idx = idx(Lia);
    
    % Prompt the user to see if there are any more bounds to set
    prompt = {'Are there other chaser metrics to constrain? 1=yes, 0=no'};
    dlg_title = 'Input';
    num_lines = 1;
    defaultans = {'0'};
    return_val = inputdlg(prompt,dlg_title,num_lines,defaultans);
    return_val = str2num(return_val{1});
end

% Plot the target space of feasible targets based on their omega and r_f
figure()
h=scatter(rad2deg(TOTAL_POINTS(7,idx)),TOTAL_POINTS(8,idx),50,TOTAL_POINTS(10,idx),'filled');
xlabel('\omega [deg/s]');
ylabel('r_f [m]')
set(gca,'XTick',OMEGA)
set(gca,'YTick',RF)
set(gca,'XTickLabel',sprintf('%3.1f\n',OMEGA))
set(gca,'YTickLabel',sprintf('%3.1f\n',RF))
xlim([min(OMEGA)*.8, max(OMEGA)*1.2]);
title({['Targets with Feasible Approaches Given']; ['User Selected Chaser Constraints']});


%% Target Properties that Enable X% Success for Given Chaser (Page 64 of lab notebook)
clc;
% Chaser Properties
vf_ub = 1;
vf_abs = abs(TOTAL_POINTS(6,:));
[~,idx_vf] = find(vf_abs < vf_ub);
amax_ub = 4;
amax_abs = abs(TOTAL_POINTS(4,:));
[~,idx_amax] = find(amax_abs < amax_ub);
deltaV_ub = 5;
deltaV_abs = abs(TOTAL_POINTS(3,:));
[~,idx_deltaV] = find(deltaV_abs < deltaV_ub);

idx_temp = intersect(idx_vf,idx_amax);
idx_chaser = intersect(idx_temp,idx_deltaV);

omega_ub = max(TOTAL_POINTS(7,idx_chaser));
omega_lb = min(TOTAL_POINTS(7,idx_chaser));
[~,idx_omega] = find((TOTAL_POINTS(7,idx_chaser) > omega_lb) & (TOTAL_POINTS(7,idx_chaser) < omega_ub));

percent_feasible = 100;%length(idx_omega)/length(idx_chaser) * 100;


delta_omega = 0.1; %rad/sec -- to be set by user
i = 1;
while percent_feasible > 95 %to be user-set
    % update bounds
    omega_ub = omega_ub + delta_omega;
    omega_lb = omega_lb - delta_omega;
    
    % find indices corresponding to the new bound
    [~,idx_omega] = find((TOTAL_POINTS(7,:) > omega_lb) & (TOTAL_POINTS(7,:) < omega_ub));
    
    % find the percent feasible with the new bound
    percent_feasible = length(idx_chaser)/length(idx_omega) * 100;
    disp([num2str(percent_feasible), '% Feasible with ', num2str(omega_lb), '<\omega [rad/s]<', num2str(omega_ub)]);
    i = i+1
end


%% Find Percent of Expected Uncertain Target Range that is Feasible with Given Chaser
clc;
% Chaser Properties
vf_ub = 0.1;
vf_abs = abs(TOTAL_POINTS(6,:));
[~,idx_vf] = find(vf_abs < vf_ub);

% Target Properties
omega = deg2rad(5);     %rotation rate [rad/s]
rf = 1;                 %docking port radius [m]

omega_sigma = deg2rad(1);   %1sigma bound on omega [rad/s]
rf_sigma = 0.01;            %1sigma bound on rf [m]

rf_ub = rf + rf_sigma;
rf_lb = rf - rf_sigma;
omega_ub = omega + omega_sigma;
omega_lb = omega - omega_sigma;

[~,idx_rf] = find((TOTAL_POINTS(8,:) > rf_lb) & (TOTAL_POINTS(8,:) < rf_ub));
[~,idx_omega] = find((TOTAL_POINTS(7,:) > omega_lb) & (TOTAL_POINTS(7,:) < omega_ub));
idx_temp = intersect(idx_rf,idx_omega);

num_targets = length(idx_temp);
num_with_valid_chasers = length(intersect(idx_temp,idx_vf));

percent_feasible = num_with_valid_chasers/num_targets * 100;
disp(['Percent of targets that the chaser limits can dock to: ', num2str(percent_feasible)]);


%% Plotting Percent of Expected Uncertain Target Range that is Feasible Over Multiple Chasers
clc;

% Target Properties
omega = deg2rad(5);     %rotation rate [rad/s]
rf = 1;                 %docking port radius [m]

omega_sigma = deg2rad(1);   %1sigma bound on omega [rad/s]
rf_sigma = 0.1;            %1sigma bound on rf [m]

rf_ub = rf + rf_sigma;
rf_lb = rf - rf_sigma;
omega_ub = omega + omega_sigma;
omega_lb = omega - omega_sigma;

[~,idx_rf] = find((TOTAL_POINTS(8,:) > rf_lb) & (TOTAL_POINTS(8,:) < rf_ub));
[~,idx_omega] = find((TOTAL_POINTS(7,:) > omega_lb) & (TOTAL_POINTS(7,:) < omega_ub));
idx_temp = intersect(idx_rf,idx_omega);

num_targets = length(idx_temp);

% Chaser Properties
vf_ub_vec = linspace(0.01,.5,100);
vf_abs = abs(TOTAL_POINTS(6,:));
for i = 1:length(vf_ub_vec)
    vf_ub = vf_ub_vec(i);
    [~,idx_vf] = find(vf_abs < vf_ub);
    num_with_valid_chasers_vf(i) = length(intersect(idx_temp,idx_vf));
    percent_feasible_vf(i) = num_with_valid_chasers_vf(i)/num_targets * 100;

end

amax_abs = abs(TOTAL_POINTS(4,:));
% amax_ub_vec = linspace(0.01,25,100);
amax_ub_vec = linspace(0.01,100,100);
for i = 1:length(amax_ub_vec)
    amax_ub = amax_ub_vec(i);
    [~,idx_amax] = find(amax_abs < amax_ub);
    num_with_valid_chasers_amax(i) = length(intersect(idx_temp,idx_amax));
    percent_feasible_amax(i) = num_with_valid_chasers_amax(i)/num_targets * 100;
end

deltaV_abs = abs(TOTAL_POINTS(3,:));
deltaV_ub_vec = linspace(0.01,5,100);
for i = 1:length(deltaV_ub_vec)
    deltaV_ub = deltaV_ub_vec(i);
    [~,idx_deltaV] = find(deltaV_abs < deltaV_ub);
    num_with_valid_chasers_deltaV(i) = length(intersect(idx_temp,idx_deltaV));
    percent_feasible_deltaV(i) = num_with_valid_chasers_deltaV(i)/num_targets * 100;
end

vf_ub_vec_normalized = vf_ub_vec./max(vf_ub_vec);
amax_ub_vec_normalized = amax_ub_vec./max(amax_ub_vec);
deltaV_ub_vec_normalized = deltaV_ub_vec./max(deltaV_ub_vec);

% This plot works, but isn't that useful, since it is only two variables
% plotted together.
figure();
[hAx] = plotyy(percent_feasible_vf,vf_ub_vec,percent_feasible_amax, amax_ub_vec);
xlabel('Feasible Percent of Targets [%]');
title('Final Velocity and Peak Acceleration Feasibility Percentage')
set(get(hAx(1), 'Ylabel'), 'String', 'Final Contact Velocity [m/s]');
set(get(hAx(2), 'Ylabel'), 'String', 'Peak Acceleration [m/s^2]');
width=550;
height=400;
set(gcf,'units','points','position',[100,100,width,height])
grid on;

figure();
plot(percent_feasible_vf,vf_ub_vec_normalized,percent_feasible_amax, amax_ub_vec_normalized, percent_feasible_deltaV,deltaV_ub_vec_normalized);
xlabel('Feasible Percent of Targets [%]');
title({['Feasibility Percentage for Normalized Chaser Capabilities']; ['v_f<', num2str(max(vf_ub_vec)), '[m/s], a_{max}<',num2str(max(amax_ub_vec)), '[m/s^2], \DeltaV<',num2str(max(deltaV_ub_vec)),'[m/s]']; [num2str(rf_lb), '<r_f [m]<', num2str(rf_ub), ', ', num2str(rad2deg(omega_lb)),'<\omega [deg/s]<',num2str(rad2deg(omega_ub))]});
legend('v_f', 'a_{max}', '\DeltaV','Location','NorthWest')
ylabel('Normalized Parameter Value')
grid on;

%% Find Feasible Chasers for a Target
% Target Properties
omega = deg2rad(5);     %rotation rate [rad/s]
rf = 1;                 %docking port radius [m]

omega_sigma = deg2rad(1);   %1sigma bound on omega [rad/s]
rf_sigma = 0.1;            %1sigma bound on rf [m]
vf_ub = .1;

idx = [];
idx_temp = [];
rf_ub = rf + rf_sigma;
rf_lb = rf - rf_sigma;
omega_ub = omega + omega_sigma;
omega_lb = omega - omega_sigma;

[~,idx_rf] = find((TOTAL_POINTS(8,:) > rf_lb) & (TOTAL_POINTS(8,:) < rf_ub));
[~,idx_omega] = find((TOTAL_POINTS(7,:) > omega_lb) & (TOTAL_POINTS(7,:) < omega_ub));
vf_abs = abs(TOTAL_POINTS(6,:));
[~,idx_vf] = find(vf_abs < vf_ub);

idx_temp = intersect(idx_rf,idx_omega);
idx = intersect(idx_temp,idx_vf);

b = TOTAL_POINTS(1,idx);
c = TOTAL_POINTS(2,idx);
dv_tot = TOTAL_POINTS(3,idx);
a_max = TOTAL_POINTS(4,idx);
tf = TOTAL_POINTS(5,idx);
vf = TOTAL_POINTS(6,idx);
omega = TOTAL_POINTS(7,idx);
rf = TOTAL_POINTS(8,idx);
r0 = TOTAL_POINTS(9,idx);
utility = TOTAL_POINTS(10,idx);

figure()
scatter3(dv_tot, a_max,r0,100,utility); hold on;
% plot3(bestPoint(3,:), bestPoint(4,:),bestPoint(9,:),'ok','MarkerSize',12);
% colormap jet
colormap winter
grid on;
title({'Utility-Colored \DeltaV vs.'; 'Peak Accleration and Starting Radius'});
xlabel('\DeltaV [m/s]');
ylabel('a_{max} [m/s^2]');
zlabel('r_0 [m]')

