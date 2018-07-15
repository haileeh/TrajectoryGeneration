%% Multi Attribut Utility Analysis - PICK BEST OPERATING POINT of R0, B, C

% Define the weights into a normalized weighting vector
% Apply negative weight as costs, positive weight as benefit
% DV, a_max, tf, vf

% Weights that must be user selected
weights = -[1 1 1 1]'; 

% Weights based on standard deviation of each metric
% weights = -[std(TOTAL_POINTS(3,~isnan(TOTAL_POINTS(3,:)))); std(TOTAL_POINTS(4,~isnan(TOTAL_POINTS(4,:)))); std(TOTAL_POINTS(5,~isnan(TOTAL_POINTS(5,:)))); std(TOTAL_POINTS(6,~isnan(TOTAL_POINTS(6,:))))]; 


% Normalize the weights to get a unit vector
weights = weights/norm(weights);

% Normalize the metrics for each point in TOTAL_POINTS
TOTAL_POINTS_NORMALIZED = TOTAL_POINTS;
TOTAL_POINTS_NORMALIZED(isnan(TOTAL_POINTS_NORMALIZED))=0;

for i = 3:6 %DV, amax, tf, vf
    TOTAL_POINTS_NORMALIZED(i,:) = TOTAL_POINTS_NORMALIZED(i,:)./norm(TOTAL_POINTS_NORMALIZED(i,:));
end
TOTAL_POINTS_NORMALIZED(6,:) = abs(TOTAL_POINTS_NORMALIZED(6,:));

% Pick the point which yields the max utility
utility = zeros(1,size(TOTAL_POINTS_NORMALIZED,2));
for i = 1:size(TOTAL_POINTS_NORMALIZED,2)
   utility(i) = sum(TOTAL_POINTS_NORMALIZED(3:6,i).*weights); 
end
utility(utility == 0) = NaN;
TOTAL_POINTS(10,:) = utility;

maxval = min(abs(TOTAL_POINTS(10,:)));
I = find(abs(TOTAL_POINTS(10,:)) == maxval);
disp(['Best Point: b=', num2str(TOTAL_POINTS(1,I)), 'm, c=', num2str(TOTAL_POINTS(2,I)), ', r0=', num2str(TOTAL_POINTS(9,I))]);
disp(['DV=', num2str(TOTAL_POINTS(3,I)), 'm/s, a_max=', num2str(TOTAL_POINTS(4,I)), 'm/s^2']);
disp(['tf=', num2str(TOTAL_POINTS(5,I)), 's, vf=', num2str(TOTAL_POINTS(6,I)),'m/s']);
disp(['utility=', num2str(TOTAL_POINTS(10,I))]);
bestPoint = TOTAL_POINTS(:,I);

%% Plot utility as a function of other variables
b = TOTAL_POINTS(1,:);
c = TOTAL_POINTS(2,:);
dv_tot = TOTAL_POINTS(3,:);
a_max = TOTAL_POINTS(4,:);
tf = TOTAL_POINTS(5,:);
vf = TOTAL_POINTS(6,:);
omega = TOTAL_POINTS(7,:);
rf = TOTAL_POINTS(8,:);
r0 = TOTAL_POINTS(9,:);
utility = TOTAL_POINTS(10,:);

figure()
scatter3(dv_tot./(omega.*rf), a_max./(omega.^2./rf),utility,9,r0); hold on;
plot3(bestPoint(3,:)./(bestPoint(7,:).*bestPoint(8,:)), bestPoint(4,:)./(bestPoint(7,:).^2.*bestPoint(8,:)),bestPoint(10,:),'ok','MarkerSize',12);
colormap jet
grid on;
title({'Trajectory Utility vs.'; 'Normalized Peak Accleration and Total \Delta V'});
xlabel('\DeltaV/(\omega r_f)');
ylabel('a_{max}/(\omega^2 r_f)');
zlabel('Utility')
%%
figure()
scatter3(dv_tot./(omega.*rf), a_max./(omega.^2./rf),r0./rf,9,utility); hold on;
% plot3(bestPoint(3,:)./(bestPoint(7,:).*bestPoint(8,:)), bestPoint(4,:)./(bestPoint(7,:).^2.*bestPoint(8,:)),bestPoint(10,:),'ok','MarkerSize',12);
colormap winter
grid on;
title({'Utility-Colored Normalized \DeltaV vs.'; 'Normalized Peak Accleration and Starting Radius'});
xlabel('\DeltaV/(\omega r_f)');
ylabel('a_{max}/(\omega^2 r_f)');
zlabel('r_0/r_f')
%%
figure()
scatter3(b, c,utility,9,utility); hold on;
plot3(bestPoint(1,:),bestPoint(2,:),bestPoint(10,:),'ok','MarkerSize',12);
colormap jet
grid on;
title('Trajectory Utility vs. Trajectory Parameters');
xlabel('b [m]');
ylabel('c');
zlabel('Utility')

figure()
scatter3(r0./rf, tf,utility,9,utility); hold on;
plot3(bestPoint(9,:)./bestPoint(8,:),bestPoint(5,:),bestPoint(10,:),'ok','MarkerSize',12);
colormap jet
grid on;
title({'Trajectory Utility vs.'; 'Normalized Starting Radius and Time to Dock'});
xlabel('r_0/r_f');
ylabel('Time to Dock [s]');
zlabel('Utility')

figure()
scatter3(vf./(omega.*rf), tf,utility,9,utility); hold on;
plot3(bestPoint(6,:)./(bestPoint(7,:).*bestPoint(8,:)),bestPoint(5,:),bestPoint(10,:),'ok','MarkerSize',12);
colormap jet
grid on;
title({'Trajectory Utility vs.'; 'Normalized Final Velocity and Time to Dock'});
xlabel('v_f/(\omega r_f)');
ylabel('Time to Dock [s]');
zlabel('Utility')

figure()
scatter3(r0, rf,utility,9,utility); hold on;
plot3(bestPoint(9,:),bestPoint(8,:),bestPoint(10,:),'ok','MarkerSize',12);
set(gcf, 'Renderer','OpenGL');
colormap jet
grid on;
title({'Trajectory Utility vs.'; 'Initial and Final Radii'});
xlabel('Starting Radius [m]');
ylabel('Final Readius [m]');
zlabel('Utility')

figure()
scatter3(b, c ,r0./rf,9,utility); hold on;
plot3(bestPoint(1,:),bestPoint(2,:),bestPoint(9,:)./bestPoint(8,:),'ok','MarkerSize',12);
colormap jet
grid on;
title({'Utility-Colored Points for Trajectory Parameters and'; 'Normalized Starting Radius'});
xlabel('b [m]');
ylabel('c');
zlabel('r_0/r_f')

figure()
scatter3(omega.*rf, vf ,dv_tot,9,utility); hold on;
plot3(bestPoint(7,:).*bestPoint(8,:),bestPoint(6,:),bestPoint(3,:),'ok','MarkerSize',12);
colormap jet
grid on;
title({'Utility-Colored Points for Docking Point Velocity,'; 'Final Radial Velocity, and \DeltaV Required'});
xlabel('\omega r_f [m/sec]');
ylabel('v_f [m/s]');
zlabel('\DeltaV [m/s]')

%% Best Point as Distance to Target Closes
clear BEST_POINTS
clear idx_best_at_r0
clear idx
for i = 1:length(r0_list)
    idx = find(TOTAL_POINTS(9,:)==r0_list(i));
    points_at_radius = TOTAL_POINTS(:,idx);
    max_utility = max(points_at_radius(10,:));
    idx_best = find(points_at_radius(10,:)==max_utility);
    idx_best_at_r0{i} = idx_best;
    TempPoints = points_at_radius(:,idx_best_at_r0{i});
    if size(TempPoints,2) > 0
        BEST_POINTS(:,i) = TempPoints(:,1);
    end
end

b = BEST_POINTS(1,:);
c = BEST_POINTS(2,:);
dv_tot = BEST_POINTS(3,:);
a_max = BEST_POINTS(4,:);
tf = BEST_POINTS(5,:);
vf = BEST_POINTS(6,:);
omega = BEST_POINTS(7,:);
rf = BEST_POINTS(8,:);
r0 = BEST_POINTS(9,:);
utility = BEST_POINTS(10,:);

figure()
subplot(3,1,[1; 2])
scatter3(b, c ,r0./rf,50,utility); hold on;
plot3(bestPoint(1,:),bestPoint(2,:),bestPoint(9,:)./bestPoint(8,:),'ok','MarkerSize',12);
colormap jet
grid on;
title({'Utility-Colored Points for Trajectory Parameters and'; 'Normalized Starting Radius'});
xlabel('b [m]');
ylabel('c');
zlabel('r_0/r_f')
subplot(3,1,3)
plot(r0,utility)
title({'Utility vs. Starting Radius'});
xlabel('r_0 [m]');
ylabel('Utility');


%% Changing Between BEST_POINTS to Follow Best Overall Path
b = BEST_POINTS(1,:);
c = BEST_POINTS(2,:);
dv_tot = BEST_POINTS(3,:);
a_max = BEST_POINTS(4,:);
tf = BEST_POINTS(5,:);
vf = BEST_POINTS(6,:);
omega = BEST_POINTS(7,:);
rf = BEST_POINTS(8,:);
r0 = BEST_POINTS(9,:);
utility = BEST_POINTS(10,:);

r0_list = r0;           %m
r_start = max(r0);      %m
r = r_start;            %m
theta = 0;              %rad
t = 0;                  %s
time_history = 0;       %s
dt = 0.001;             %s
i = 1;
b_now = b(end);         %m
c_now = c(end);
omega = mean(omega);    %rad/s

while r > min(rf)
    i = i+1;
    t(i) = t(i-1) + dt;
    time_history(i) = time_history(i-1) + dt;
    theta(i) = theta(i-1) + omega*dt;
    r(i) = r_start - b_now(i-1)*(omega*t(i))^(1/c_now(i-1));
    if r(i) < r0_list(end)
        if (length(r0_list) > 1) && (length(b_now) > 1)
            idx = length(r0_list);
            b_now(i) = b(idx-1);
            c_now(i) = c(idx-1);
            r_start = r(i);
            r0_list(end) = []
            t(i) = 0;
        else
            b_now(i) = b_now(i-1);
            c_now(i) = c_now(i-1);
        end
    else
        b_now(i) = b_now(i-1);
        c_now(i) = c_now(i-1);
    end
end


figure();
subplot(3,1,[1,2])
polar(theta,r);
% ylim([-1, 10])
% xlim([-1, 10])
title('Archimedes Spiral with Changing Best b,c');
subplot(3,1,3)
plot(time_history,b_now,'b'); hold on;
plot(time_history,c_now,'--r');
plot(time_history,r,'-.k');
title('Trajectory Parameters Over Time')
legend('b','c','r')
xlabel('Time [s]')
ylabel('Value; b,r in [m]')

% figure()
% plot(t,theta)

