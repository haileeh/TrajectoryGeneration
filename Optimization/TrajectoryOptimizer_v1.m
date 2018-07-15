%% Set-up
%clear all; 
clc;
close all; 

%% Initialize cells and struct arrays
fit_error_exponential = [];
computational_time_exponential = [];
computational_DV_exponential = [];
run_num = 1;
run_num_exponential = 1;
RESULTS_exponential = [];
optimal_time_history = [];
optimal_radial_distribution = [];

%% Set variables
num_dt = 100;
NUM_TERMS = 6;
omega = 5/norm([0 0 1])*[0 0 1]';%5/norm([0 0.25 1])*[0 0.25 1]'; %5/norm([0 0 1])*[0 0 1]';
r0 = 10;
rf = 1;
inertia_ratio = [1 1 1]';%[2 2 1]';
fit_type = {'exponential'}; %{'polynomial','exponential','logistic'}
%Tend = 30; %120;%31.3;%180; %60;
%Tend = [1; 20; 33; 55; 240];
Tend = 33;%[30 33 40];
for j=1:length(Tend)
%% 1. Find the Fuel Optimal Trajectory

% Setup for fmincon call
A = ones(num_dt);
b = Tend(j)*ones(num_dt,1);
Aeq = ones(1,num_dt); %[];
beq = Tend(j); %[]
lb = zeros(num_dt,1);
ub = Tend(j)*ones(num_dt,1);

% Initial guess = uniform time step
x0 = Tend(j)/num_dt * ones(num_dt,1);
%changed max's to 100,000 from 1mil
options.MaxFunEvals = 10000;
options.TolFun = 1e-15;
options.Display = 'off';
options = optimset('TolFun',1e-8,'MaxFunEvals',1000000,'MaxIter',10000000,'TolX',1e-8,'TolCon',1e-8,'Algorithm','sqp');
%options = optimoptions(@fmincon,'Algorithm','sqp','ConstraintTolerance',1e-8,...
%    'MaxFunctionEvaluations',1e5,'MaxIterations',1e5,'FunctionTolerance',...
%    1e-8, 'StepTolerance',1e-8);

f = @(dt_list)DeltaV_dt(dt_list,omega,r0,rf,inertia_ratio);
% only 5 arguments means log spacing!
[dt_out, DV_out,exitflag,output] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,[],options);
disp(['Exitflag: ',num2str(exitflag)]);
disp(['Minimum DV=',num2str(DV_out),'m/s'])

% -----Optimal Time History-----
t_out = zeros(length(dt_out),1);
for i = 1:length(dt_out)
    t_out(i+1) = t_out(i) + dt_out(i);
end

r = logspace(log10(r0),log10(rf),length(t_out)); r(end) = rf;

figure();
plot(t_out,r,'-*b');
title({'Optimal Radial Distance Profile';['\DeltaV=',num2str(DV_out),'m/s']})
xlabel('Time [s]')
ylabel('Radius [m]')

% if t_out(end) > (t_out(end-1) + 5)
%     t_out = t_out(1:end-1);
%     dt_out = dt_out(1:end-1);
% end

optimal_time_history(:,i) = t_out; x = t_out;
optimal_radial_distribution(i,:) = r;

DV_computed = DeltaV_dt_r_omega0(dt_out,optimal_radial_distribution(i,:),omega,inertia_ratio,1)
% Create struct
DV_computed_struct.dt_out(:,j) = dt_out;
DV_computed_struct.DV(j) = DV_computed;
DV_computed_struct.DV_out(j) = DV_out;
DV_computed_struct.Tend(j) = Tend(j);
DV_computed_struct.Alin(:,:,j) = Alin;
DV_computed_struct.Acen(:,:,j) = Acen;
DV_computed_struct.Acor(:,:,j) = Acor;
DV_computed_struct.Aang(:,:,j) = Aang;
end

keyboard
%% 2. Fit a parameterized approximation to the fuel optimal trajectory
[fcnfit_opt,fit_error,B_list]=fit_param_approx(optimal_radial_distribution, optimal_time_history, r0, rf, fit_type,NUM_TERMS, num_dt);
DV_fit1 = DeltaV_dt_r_omega0(dt_out,fcnfit_opt(2,:),omega,inertia_ratio,1)
keyboard
%% 3. Find the 2-term Fuel Optimal Exponential Trajectory

t_vec = linspace(0,Tend,num_dt+1);
dt = diff(t_vec);

% Nonlinear constrainted optimization function with argument b - the coefficients.
nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(t_vec,b,r0,rf,'exponential');

scaling_vec = [1; 1; 1; 1]; % does this actually do anything

% f is the objective function!
f = @(exp_coefficients)DeltaV_dt(dt,omega,r0,rf,inertia_ratio,exp_coefficients./scaling_vec);
% x0 = [8.37882458523862;0.0551170487169159;1.62117541476138;0.00268646208704073];
% why this one below selected?
%x0 = [19.9995004517454;0.0487356457859513;-9.99950045174537;0.0337307942865192];
% x0 = [0 0 5 0]';
% x0 = [9; 0; 2; 0];
x0 = B_list{2}; %selecting results from step 2 as initial guess

opts = optimoptions(@fmincon,'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',f,'x0',x0,'Aeq',...
    repmat([1,0],1,2),'beq',r0,'lb',[],'ub',[],'nonlcon',nonlcon,'options',opts);
gs = GlobalSearch;
[exp_coeff_out,DV_out2] = run(gs,problem);
% exp_coeff_out is the value yielding DV_out (global minimum)

disp(['Minimum DV=',num2str(DV_out2),'m/s'])

t_out = zeros(length(dt),1);
for i = 1:length(dt)
    t_out(i+1) = t_out(i) + dt(i);
end

% Calculates radial profile of Chaser by time history done for optimization
% of the coefficients - 2-term exponential
r = exp_coeff_out(1)*exp(-exp_coeff_out(2)*t_out) + ...
    exp_coeff_out(3)*exp(-exp_coeff_out(4)*t_out);

x = t_out;  %Re-naming variable
yx = r;     %Re-naming variable; this is the optimized exp trajectory
% This is the 100 points along optimal 2-term exponential trajectory
figure();
plot(t_out,r);
title({'Optimal Radial Distance Profile';['\DeltaV=',num2str(DV_out2),'m/s']})
xlabel('Time [s]')
ylabel('Radius [m]')

% This is finding the deltaV from the optimal exp trajectory
DeltaV_Optimal = DeltaV_dt_r_omega0(dt,yx',omega,inertia_ratio,1);

keyboard
%% 4. Fitting the exponential optimal trajectory!
[fcnfit_exp_opt,fit_error,B_list]=fit_param_approx(yx', x, r0, rf, fit_type,NUM_TERMS, num_dt);
% perhaps add option for specifiying which optimal trajecotry this is
num_terms = NUM_TERMS;
for i = 1:num_terms
    switch i
        case 1
            fit_coeffs = B_list{1}(:,end);%B1(:,end);
            r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x);
        case 2
            fit_coeffs =  B_list{2}(:,end);%B2(:,end);
            r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x);
        case 3
            fit_coeffs =  B_list{3}(:,end);%B3(:,end);
            r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x) + fit_coeffs(5).*exp(-fit_coeffs(6).*x);
        case 4
            fit_coeffs =  B_list{4}(:,end);%B4(:,end);
            r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x) + fit_coeffs(5).*exp(-fit_coeffs(6).*x) +...
                fit_coeffs(7).*exp(-fit_coeffs(8).*x);
        case 5
            fit_coeffs =  B_list{5}(:,end);%B5(:,end);
            r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x) + fit_coeffs(5).*exp(-fit_coeffs(6).*x) +...
                fit_coeffs(7).*exp(-fit_coeffs(8).*x) + fit_coeffs(9).*exp(-fit_coeffs(10).*x);
        case 6
            fit_coeffs =  B_list{6}(:,end);%B6(:,end);
            r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x) + fit_coeffs(5).*exp(-fit_coeffs(6).*x) +...
                fit_coeffs(7).*exp(-fit_coeffs(8).*x) + fit_coeffs(9).*exp(-fit_coeffs(10).*x) + fit_coeffs(11).*exp(-fit_coeffs(12).*x);
    end
    
    R_FIT2(:,i) = r_fit; %fcnfit(i,:);
    DeltaV_Fit(i) = DeltaV_dt_r_omega0(dt,r_fit',omega,inertia_ratio,0);%spiralDeltaV_dt_r_omega0(dt,fcnfit(i,:),omega,inertia_ratio,0);
    
end

DeltaV_Fit_Normalized = DeltaV_Fit/DeltaV_Optimal; % optimal exp trajectory

% where did this come from? I think David just decided the two-term fit was
% best, so he implemented it here. are these coefs from B2?
% fit_coeffs = [8.37882458523862;0.0551170487169159;1.62117541476138;0.00268646208704073];
% r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x);
% DeltaV_TwoTermBestFit = DeltaV_dt_r_omega0(dt,r_fit',omega,inertia_ratio,0)

figure()
plot(x,R_FIT2(:,2)); % but R_FIT has 6 columns? so what is actually being plotted?
grid on;
xlabel('Time [s]')
ylabel('Radius [m]')
title('Radial Distance of Exponential Fits')

figure()
plot(DeltaV_Fit_Normalized);
grid on;
xlabel('Number of Terms')
ylabel('Normalized \DeltaV')
title('Fit \DeltaV Normalized by Optimal \DeltaV')

DeltaV_dt_r_omega0(dt,R_FIT2(:,2)',omega,inertia_ratio,1)
keyboard

%% Plot all 4
figure
plot(optimal_time_history,optimal_radial_distribution,'-*b'); hold on;
plot(optimal_time_history,fcnfit_opt(2,:),'-ok');
plot(t_out, yx,'-.g');
plot(t_out, fcnfit_exp_opt(2,:),'--r'); grid on;
xlabel('Time [s]');
ylabel('Radius [m]');
legend('100-Point Optimal Trajectory','2-Term Fit of Optimal','100 Points Along Optimal 2-Term Exponential Trajectory',...
    '2-Term Fit to 2-Term Exponential Trajectory');

%% Compare with diffeq
% make the propagator into a function?