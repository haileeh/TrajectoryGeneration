clear all; clc;
close all;

R_LIST_polynomial = {};
T_LIST_polynomial = {};
DV_LIST_polynomial = [];
OMEGA_LIST_polynomial = [];
fit_error_polynomial = [];
R_LIST_exponential = {};
T_LIST_exponential = {};
DV_LIST_exponential = [];
OMEGA_LIST_exponential = [];
fit_error_exponential = [];
computational_time_exponential = [];
computational_DV_exponential = [];

num_dt_list = 100
NUM_TERMS = 6;
num_terms_best_exponential = 2;
% Tend_list = [1 2 3 4 5] * 60
% omega_list = [0 0 5;...
%     0 0 1;...
%     0 5 0;...
%     0 1 0;...
%     5 0 0;...
%     1 0 0;...
%     5 5 5;...
%     1 1 1]'
% omega_list = [0 0 1;...
%     0 0 2;...
%     0 0 3;...
%     0 0 4;...
%     0 0 5;...
%     0 0 6;...
%     0 0 7;...
%     0 0 8;...
%     0 0 9;...
%     0 0 10]'
% omega_list = [0 0 5;...
%     0 0 6;...
%     0 0 7;...
%     0 0 8]'
omega_list = 5/norm([0 0 1])*[0 0 1]'
% omega_list = 5/norm([0 0.25 1])*[0 0.25 1]'
% omega_list = 5/norm([1 1 1])*[1 1 1]'
% % omega_list = [5/norm([0 0 1])*[0 0 1]',...
% %     5/norm([1 1 1])*[1 1 1]',...
% %     5/norm([0 0.05 1])*[0 0.05 1]'];  
% r0_list = [2 4 5 6 7 8 9 10] %m
% rf_list = [0.5 0.6 0.7 0.8 0.9 1] %m
% inertia_ratio_list = [1 1 1;...
%     1.2 1 1;...
%     1 1.2 1;...
%     1 1 1.2;...
%     1.2 1.2 1;...
%     1 1.2 1.2;...
%     1.2 1 1.2;...
%     1.5 1.2 1;...
%     1.2 1.5 1;...
%     1 1.2 1.5;...
%     1 1.5 1.2;...
%     1.5 1 1.2;...
%     1.2 1 1.5]'
fit_type_list = {'exponential'}%{'polynomial','exponential','logistic'}
Tend_list = [1 20 33 55 240]; %60%[30 60 120]
% omega_list = [0 0 1;...
%     0 0 5]'
r0_list = 10%[5 10 20]
rf_list = 1%[.5 0.75 1]
% % inertia_ratio_list = [1 1 1;...
% %     1.2 1 1;...
% %     2 2 1]'
% inertia_ratio_list = [1 1 1;...
%     1.2 1 1]'
inertia_ratio_list = [1 1 1]'
% inertia_ratio_list = [2 2 1]'
% inertia_ratio_list = [2 0.5 1]'
% inertia_ratio_list = [.35 .35 1]'
% inertia_ratio_list = [2 2 1]'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE SURE PLOT TWO CELLS BELOW HAS SAME LIST CALLED OUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TotalRuns = length(num_dt_list) * length(Tend_list) * size(omega_list,2) * length(r0_list) * length(rf_list) * size(inertia_ratio_list,2) * length(fit_type_list)

run_num = 1;
run_num_exponential = 1;
run_num_polynomial = 1;
RESULTS_polynomial = [];
RESULTS_exponential = [];
optimal_time_history = [];
optimal_radial_distribution = [];

% Compute the Best Fit for Each Trajectory
for num_dt_idx = 1:length(num_dt_list)
    for Tend_idx = 1:length(Tend_list)
        for omega_idx = 1:size(omega_list,2)
            for r0_idx = 1:length(r0_list)
                for rf_idx = 1:length(rf_list)
                    for inertia_idx = 1:size(inertia_ratio_list,2)
                        
                        tic
                        
                        inertia_ratios = inertia_ratio_list(:,inertia_idx);
                        rf = rf_list(rf_idx)
                        r0 = r0_list(r0_idx)
                        omega = omega_list(:,omega_idx)
                        Tend = Tend_list(Tend_idx)
                        num_dt = num_dt_list(num_dt_idx);
                        fcnfit = zeros(NUM_TERMS,num_dt+1);
                        if TotalRuns == 1
                            plot_single = 1;
                        else
                            plot_single = 0;
                        end
                        A = ones(num_dt);
                        b = Tend*ones(num_dt,1);
                        Aeq = ones(1,num_dt);
                        beq = Tend;
                        lb = zeros(num_dt,1);
                        ub = Tend*ones(num_dt,1);
                        
                        x0 = Tend/num_dt * ones(num_dt,1);
                        options.MaxFunEvals = 10000;
                        options.TolFun = 1e-15;
                        options.Display = 'off';
                        options = optimset('TolFun',1e-8,'MaxFunEvals',1000000,'MaxIter',10000000,'TolX',1e-8,'TolCon',1e-8,'Algorithm','sqp');
                        
                        f = @(dt_list)spiralDeltaV_dt_VariableStudy(dt_list,omega,r0,rf,inertia_ratios);
                        [dt_out, DV_out,exitflag,output] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,[],options);
                        disp(['Exitflag: ',num2str(exitflag)]);
                        disp(['Minimum DV=',num2str(DV_out),'m/s'])
                        
                        exp_coeff_out = zeros(4,1);
                        t_vec = linspace(0,Tend,num_dt+1);%(logspace(log10(0.0001),log10(Tend),num_dt+1));
                        dt_list = diff(t_vec);
                        nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(t_vec,b,r0,rf,'exponential');
                        if num_terms_best_exponential == 2
                            scaling_vec = [1; 1; 1; 1];
                            f = @(exp_coefficients)spiralDeltaV_dt_VariableStudy(dt_list,omega,r0,rf,inertia_ratios,exp_coefficients./scaling_vec);
%                             x0 = [8.37882458523862;0.0551170487169159;1.62117541476138;0.00268646208704073];
                            x0 = [19.9995004517454;0.0487356457859513;-9.99950045174537;0.0337307942865192];
%                           x0 = [0 0 5 0]';
%                           x0 = [9; 0; 2; 0];
%                           [exp_coeff_out, DV_out,exitflag,output] = fmincon(f,x0,[],[],repmat([1,0],1,2),r0,[],[],nonlcon,options);
                        
                            opts = optimoptions(@fmincon,'Algorithm','interior-point');
                            problem = createOptimProblem('fmincon','objective',f,'x0',x0,'Aeq',repmat([1,0],1,2),'beq',r0,'lb',[],'ub',[],'nonlcon',nonlcon,'options',opts);
                            gs = GlobalSearch;
                            [exp_coeff_out,DV_out] = run(gs,problem);
                        elseif num_terms_best_exponential == 3
                            scaling_vec = [1; 1; 1; 1; 1; 1];
                            f = @(exp_coefficients)spiralDeltaV_dt_VariableStudy(dt_list,omega,r0,rf,inertia_ratios,exp_coefficients./scaling_vec);
                            x0 = [9.04079194195414;0.0509358869930513;2.16063608354774;-0.00805623257712262;-1.20142802550188;-0.0106790005587400];
%                           [exp_coeff_out, DV_out,exitflag,output] = fmincon(f,x0,[],[],repmat([1,0],1,2),r0,[],[],nonlcon,options);
                        
                            opts = optimoptions(@fmincon,'Algorithm','interior-point');
                            problem = createOptimProblem('fmincon','objective',f,'x0',x0,'Aeq',repmat([1,0],1,3),'beq',r0,'lb',-10*ones(6,1),'ub',10*ones(6,1),'nonlcon',nonlcon,'options',opts);
                            gs = GlobalSearch;
                            [exp_coeff_out,DV_out] = run(gs,problem);
                            exp_coeff_out = x0;
                        end

                        dt_out = dt_list;
%                         disp(['Exitflag: ',num2str(exitflag)]);
                        disp(['Minimum DV=',num2str(DV_out),'m/s'])
                        
                        
                        
                        t_out = 0;
                        for i = 1:length(dt_out);
                            t_out(i+1) = t_out(i) + dt_out(i);
                        end
                        
                        if sum(exp_coeff_out) == 0
                            r = logspace(log10(10),log10(1),length(t_out)); r(end) = rf;
                        elseif num_terms_best_exponential == 2
                            r = exp_coeff_out(1)*exp(-exp_coeff_out(2)*t_out) + ...
                            exp_coeff_out(3)*exp(-exp_coeff_out(4)*t_out);
                        elseif num_terms_best_exponential == 3
                            r = exp_coeff_out(1)*exp(-exp_coeff_out(2).*t_out) + ...
                            exp_coeff_out(3)*exp(-exp_coeff_out(4).*t_out) + ...
                            exp_coeff_out(5)*exp(-exp_coeff_out(6).*t_out);
                        end
                        
                        x = t_out;
%                         x = smooth(linspace(1,length(t_out),length(t_out)),t_out,0.1,'rloess')'; x(1) = 0; x(end) = t_out(end);
%                         dt_out = smooth(linspace(1,length(dt_out),length(dt_out)),dt_out,0.1,'rloess')';
                        yx = r;
                        optimal_time_history = [optimal_time_history; x];
                        optimal_radial_distribution = [optimal_radial_distribution; yx];
                        
                        if plot_single == 1
                            figure();
                            plot(t_out);
                            title('Time over Index');
                            xlabel('Index')
                            ylabel('Time [s]');
                            
                            figure();
                            plot(t_out,r);
                            title({'Optimal Radial Distance Profile';['\DeltaV=',num2str(DV_out),'m/s']})
                            xlabel('Time [s]')
                            ylabel('Radius [m]')
                        end
                        
                        for fit_type_idx = 1:length(fit_type_list)
                            % Fit an exponential with as many terms as necessary
                            % https://www.mathworks.com/matlabcentral/answers/124680-3-term-exponential-fit
                            B1 = []; B2 = []; B3 = []; B4 = []; B5 = []; B6 = []; B7 = []; B8 = []; B9 = []; B10 = []; B11 = []; B12 = [];
                            for num_terms = 1:NUM_TERMS
                                %                         fit_type = 'exponential';
                                %                         fit_type = 'polynomial';
                                %                         fit_type = 'exponential_with_r0';
                                fit_type = fit_type_list{fit_type_idx};
                                if strcmp(fit_type,'exponential')
                                    switch num_terms
                                        case 1
                                            y = @(b,x) b(1).*exp(-b(2).*x);
                                        case 2
                                            y = @(b,x) b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x);
                                        case 3
                                            y = @(b,x) b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x) + b(5).*exp(-b(6).*x);
                                        case 4
                                            y = @(b,x) b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x) + b(5).*exp(-b(6).*x) + b(7).*exp(-b(8).*x);
                                        case 5
                                            y = @(b,x) b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x) + b(5).*exp(-b(6).*x) + b(7).*exp(-b(8).*x) + b(9).*exp(-b(10).*x);
                                        case 6
                                            y = @(b,x) b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x) + b(5).*exp(-b(6).*x) + b(7).*exp(-b(8).*x) + b(9).*exp(-b(10).*x) + b(11).*exp(-b(12).*x);
                                        otherwise
                                            disp('Need to create new exponential function in switch statement');
                                    end
                                elseif strcmp(fit_type,'logistic')
                                    %b(1)=C  b(2)=v  b(3)=Q  b(4)=b  b(5)=A  b(6)=k
                                    %y = A + (k-A)./(C+Q*exp(-b.*x)).^(1/v)
                                    y = @(b,x) b(5) + (b(6)-b(5))./(b(1)+b(3)*exp(-b(4).*x)).^(1/b(2));
                                elseif strcmp(fit_type,'exponential_with_r0')
                                    switch num_terms
                                        case 1
                                            y = @(b,x) r(1) - (b(1).*exp(b(2).*x));
                                        case 2
                                            y = @(b,x) r(1) - (b(1).*exp(b(2).*x) + b(3).*exp(b(4).*x));
                                        case 3
                                            y = @(b,x) r(1) - (b(1).*exp(b(2).*x) + b(3).*exp(b(4).*x) + b(5).*exp(b(6).*x));
                                        case 4
                                            y = @(b,x) r(1) - (b(1).*exp(b(2).*x) + b(3).*exp(b(4).*x) + b(5).*exp(b(6).*x) + b(7).*exp(b(8).*x));
                                        case 5
                                            y = @(b,x) r(1) - (b(1).*exp(b(2).*x) + b(3).*exp(b(4).*x) + b(5).*exp(b(6).*x) + b(7).*exp(b(8).*x) + b(9).*exp(b(10).*x));
                                        case 6
                                            y = @(b,x) r(1) - (b(1).*exp(b(2).*x) + b(3).*exp(b(4).*x) + b(5).*exp(b(6).*x) + b(7).*exp(b(8).*x) + b(9).*exp(b(10).*x) + b(11).*exp(b(12).*x));
                                        otherwise
                                            disp('Need to create new exponential function in switch statement');
                                    end
                                elseif strcmp(fit_type,'polynomial')
                                    switch num_terms
                                        case 1
                                            y = @(b,x) b(1).*x.^0;
                                        case 2
                                            y = @(b,x) b(2).*x.^0 + b(1).*x.^1;
                                        case 3
                                            y = @(b,x) b(3).*x.^0 + b(2).*x.^1 + b(1).*x.^2;
                                        case 4
                                            y = @(b,x) b(4).*x.^0 + b(3).*x.^1 + b(2).*x.^2 + b(1).*x.^3;
                                        case 5
                                            y = @(b,x) b(5).*x.^0 + b(4).*x.^1 + b(3).*x.^2 + b(2).*x.^3 + b(1).*x.^4;
                                        case 6
                                            y = @(b,x) b(6).*x.^0 + b(5).*x.^1 + b(4).*x.^2 + b(3).*x.^3 + b(2).*x.^4 + b(1).*x.^5;
                                        case 7
                                            y = @(b,x) b(7).*x.^0 + b(6).*x.^1 + b(5).*x.^2 + b(4).*x.^3 + b(3).*x.^4 + b(1).*x.^5 + b(1).*x.^6;
                                        case 8
                                            y = @(b,x) b(8).*x.^0 + b(7).*x.^1 + b(6).*x.^2 + b(5).*x.^3 + b(4).*x.^4 + b(3).*x.^5 + b(2).*x.^6 + b(1).*x.^7;
                                        case 9
                                            y = @(b,x) b(9).*x.^0 + b(8).*x.^1 + b(7).*x.^2 + b(6).*x.^3 + b(5).*x.^4 + b(4).*x.^5 + b(3).*x.^6 + b(2).*x.^7 + b(1).*x.^8;
                                        case 10
                                            y = @(b,x) b(10).*x.^0 + b(9).*x.^1 + b(8).*x.^2 + b(7).*x.^3 + b(6).*x.^4 + b(5).*x.^5 + b(4).*x.^6 + b(3).*x.^7 + b(2).*x.^8 + b(1).*x.^9;
                                        case 11
                                            y = @(b,x) b(11).*x.^0 + b(10).*x.^1 + b(9).*x.^2 + b(8).*x.^3 + b(7).*x.^4 + b(6).*x.^5 + b(5).*x.^6 + b(4).*x.^7 + b(3).*x.^8 + b(2).*x.^9 + b(1).*x.^10;
                                        case 12
                                            y = @(b,x) b(12).*x.^0 + b(11).*x.^1 + b(10).*x.^2 + b(9).*x.^3 + b(8).*x.^4 + b(7).*x.^5 + b(6).*x.^6 + b(5).*x.^7 + b(4).*x.^8 + b(3).*x.^9 + b(2).*x.^10 + b(1).*x.^11;
                                        otherwise
                                            disp('Need to create new polynomial function in switch statement');
                                    end %end of switch over num_terms
                                end
                                
                                
                                OLS = @(b) sum((y(b,x) - yx).^2);          % Ordinary Least Squares cost function
                                
                                for num_evals_iters = [1E10, 1E100] %num_evals_iters = [1E10, 1E50, 1E100, 1E200, 1E300, 1E400, 1E500, 1E600, 1E700, 1E800, 1E900]
                                    opts = optimset('MaxFunEvals',num_evals_iters, 'MaxIter',num_evals_iters,'TolFun',1e-7,'TolCon',1e-9,'Display', 'off');
                                    %                             fit_method = 'use_fminsearch';
                                    %                             fit_method = 'use_fmincon';
                                    if strcmp(fit_type,'polynomial')
%                                         fit_method = 'use_polyfit';
                                        fit_method = 'use_fmincon';
                                    else
                                        fit_method = 'use_fmincon';
                                    end
                                    if strcmp(fit_method,'use_fimsearch')
                                        switch num_terms
                                            case 1
                                                B = fminsearch(OLS, rand(2,1), opts);
                                                B1 = [B1,B];
                                            case 2
                                                B = fminsearch(OLS, [B_list{1};0;0], opts);
                                                B2 = [B2,B];
                                            case 3
                                                B = fminsearch(OLS, [B_list{2};0;0], opts);
                                                B3 = [B3,B];
                                            case 4
                                                B = fminsearch(OLS, [B_list{3};0;0], opts);
                                                B4 = [B4,B];
                                            case 5
                                                B = fminsearch(OLS, [B_list{4};0;0], opts);
                                                B5 = [B5,B];
                                            case 6
                                                B = fminsearch(OLS, [B_list{5};0;0], opts);
                                                B6 = [B6,B];
                                            otherwise
                                                disp('Need to create new exponential function in switch statement');
                                        end
                                    elseif strcmp(fit_method,'use_fmincon')
                                        if strcmp(fit_type,'logistic')
                                            A_nu = zeros(6,6); A_nu(2,2) = 1;
                                            Aeq_logistic = zeros(6,6); Aeq_logistic(5,5)=1; Aeq_logistic(6,6)=1;
                                            beq_logistic = [0;0;0;0;10;1];
                                            nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_GLF(x,b);
                                            B = fmincon(OLS,[1;0.1;1;0.1;10;1],[],[],[],[],[-inf, 0, -inf, -inf, 9, 0]',[inf, inf, inf, inf, inf, 2]',nonlcon,opts);
%                                             B = fmincon(OLS,[1;0.1;1;0.1;10;1],[],[],[],[],[],[],nonlcon,opts);
%                                             B=[1;0.1;1;0.1;10;1];
                                            B1 = [B1,B];
                                        elseif strcmp(fit_type,'exponential')
                                            switch num_terms
                                                case 1
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, rand(2,1),[], [],repmat([1,0],1,num_terms),r0,-10*ones(2,1),[],nonlcon,opts);
                                                    B1 = [B1,B];
                                                case 2
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{1};0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(4,1),[],nonlcon,opts);
                                                    B2 = [B2,B];
                                                case 3
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{2};0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(6,1),[],nonlcon,opts);
                                                    B3 = [B3,B];
                                                case 4
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{3};0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(8,1),[],nonlcon,opts);
                                                    B4 = [B4,B];
                                                case 5
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{4};0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(10,1),[],nonlcon,opts);
                                                    B5 = [B5,B];
                                                case 6
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{5};0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(12,1),[],nonlcon,opts);
                                                    %                                                 B = fmincon(OLS, [B_list{5};0;0],[], [],ones(1,12),r0,-10*ones(12,1),[],[],opts);
                                                    B6 = [B6,B];
                                            end
                                        else
                                            switch num_terms
                                                case 1
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, rand(1,1),[], [],[],[],[],[],nonlcon,opts);
                                                    B1 = [B1,B];
                                                case 2
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{1};0],[], [],[],[],[],[],nonlcon,opts);
                                                    B2 = [B2,B];
                                                case 3
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{2}';0],[], [],[],[],[],[],nonlcon,opts);
                                                    B3 = [B3,B];
                                                case 4
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{3}';0],[], [],[],[],[],[],nonlcon,opts);
                                                    B4 = [B4,B];
                                                case 5
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{4}';0],[], [],[],[],[],[],nonlcon,opts);
                                                    B5 = [B5,B];
                                                case 6
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{5}';0],[], [],[],[],[],[],nonlcon,opts);
                                                    B6 = [B6,B];
                                                case 7
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{6}';0],[], [],[],[],[],[],nonlcon,opts);
                                                    B7 = [B7,B];
                                                case 8
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{7}';0],[], [],[],[],[],[],nonlcon,opts);
                                                    B8 = [B8,B];
                                                case 9
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{8}';0],[], [],[],[],[],[],nonlcon,opts);
                                                    B9 = [B9,B];
                                                case 10
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{9}';0],[], [],[],[],[],[],nonlcon,opts);
                                                    B10 = [B10,B];
                                                case 11
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{10}';0],[], [],[],[],[],[],nonlcon,opts);
                                                    B11 = [B11,B];
                                                case 12
                                                    %                                                 nonlcon = @rf_con;
                                                    nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                                                    B = fmincon(OLS, [B_list{11}';0],[], [],[],[],[],[],nonlcon,opts);
                                                    B12 = [B12,B];
                                            end
                                        end
                                    elseif strcmp(fit_method,'use_polyfit')
                                        switch num_terms
                                            case 1
                                                B = polyfit(x,yx,0);
                                                B1 = [B1,B];
                                            case 2
                                                B = polyfit(x,yx,1);
                                                B2 = [B2,B];
                                            case 3
                                                B = polyfit(x,yx,2);
                                                B3 = [B3,B];
                                            case 4
                                                B = polyfit(x,yx,3);
                                                B4 = [B4,B];
                                            case 5
                                                B = polyfit(x,yx,4);
                                                B5 = [B5,B];
                                            case 6
                                                B = polyfit(x,yx,5);
                                                B6 = [B6,B];
                                            case 7
                                                B = polyfit(x,yx,6);
                                                B7 = [B7,B];
                                            case 8
                                                B = polyfit(x,yx,7);
                                                B8 = [B8,B];
                                            case 9
                                                B = polyfit(x,yx,8);
                                                B9 = [B9,B];
                                            case 10
                                                B = polyfit(x,yx,9);
                                                B10 = [B10,B];
                                            case 11
                                                B = polyfit(x,yx,10);
                                                B11 = [B11,B];
                                            case 12
                                                B = polyfit(x,yx,11);
                                                B12 = [B12,B];
                                        end
                                    end
                                    if sum(fcnfit(num_terms,:)) == 0
                                        fcnfit(num_terms,:) = y(B,x);                   % Calculate function with estimated parameters
                                        fit_error(num_terms,:) = [num_terms, norm(abs(fcnfit(num_terms,:)-yx)), num_evals_iters];
                                    end
                                    error_in_this_fit = norm(abs(fcnfit(num_terms,:)-yx));
                                    if error_in_this_fit < fit_error(num_terms,2)
                                        fcnfit(num_terms,:) = y(B,x);
                                        fit_error(num_terms,:) = [num_terms, error_in_this_fit, num_evals_iters];
                                    end
                                    
                                    if strcmp(fit_type,'polynomial')
                                        RESULTS_polynomial{run_num_polynomial,num_terms} = [r0;rf;omega;Tend;error_in_this_fit;reshape(B,length(B),1);inertia_ratios];
                                        B_list{num_terms} = B';
                                        fit_error_polynomial = [fit_error_polynomial;fit_error(num_terms,:)];
                                        disp('in poly')
                                    else
                                        RESULTS_exponential{run_num_exponential,num_terms} = [r0;rf;omega;Tend;error_in_this_fit;B;inertia_ratios];
                                        B_list{num_terms} = B;
                                        fit_error_exponential = [fit_error_exponential;fit_error(num_terms,:)];
                                        disp('in exp')
                                    end
                                                                      
                                    
                                end %loop over num_evals_iters
                                if plot_single == 1
                                    figure()
                                    plot(x, yx, '*-b')
                                    hold on
                                    plot(x, fcnfit(num_terms,:), '-r')
                                    hold off
                                    grid
                                    title({'Radial Distance Exponential Fit with';[num2str(num_terms),' Terms; Norm(Error)=',num2str(norm(fcnfit(num_terms,:)-yx))]})
                                    xlabel('Time [s]')
                                    ylabel('Radial Distance [m]')
                                end
                                
                            end %loop over num_terms
                            if plot_single == 1
                                figure();
                                plot(fit_error(:,1),fit_error(:,2))
                                title('Fit Error over Number of Exponential Terms')
                                xlabel('Number of Terms')
                                ylabel('Norm of Error')
                            end
                            
                            if strcmp(fit_type,'polynomial')
                                B_LIST_polynomial{run_num_polynomial} = B_list;
                                R_LIST_polynomial{end+1} = r;
                                T_LIST_polynomial{end+1} = t_out;
                                DV_LIST_polynomial(end+1) = DV_out;
                                OMEGA_LIST_polynomial(:,end+1) = omega;
                                run_num_polynomial = run_num_polynomial + 1;
                            else
                                B_LIST_exponential{run_num_exponential} = B_list;
                                R_LIST_exponential{end+1} = r;
                                T_LIST_exponential{end+1} = t_out;
                                DV_LIST_exponential(end+1) = DV_out;
                                OMEGA_LIST_exponential(:,end+1) = omega;
                                run_num_exponential = run_num_exponential + 1;
                            end
                            disp(['Percent Complete: ',num2str(100*run_num/TotalRuns)]);
                            run_num = run_num + 1;
                            computational_time_exponential = [computational_time_exponential; toc];
                            DV_computed = spiralDeltaV_dt_r_omega0(dt_out,yx,omega_list,inertia_ratios,0)
                            computational_DV_exponential = [computational_DV_exponential; DV_computed];
                            
                        end %loop over fit_type_idx
                    end %loop over inertia ratios
                end %loop over rf
            end %loop over r0
        end %loop over omega
    end %loop over Tend
end %loop over num_dt_list

if strcmp(fit_type,'exponential')
    fit_coeffs = B6(:,end);
    r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x) + fit_coeffs(5).*exp(-fit_coeffs(6).*x) +...
        fit_coeffs(7).*exp(-fit_coeffs(8).*x) + fit_coeffs(9).*exp(-fit_coeffs(10).*x) + fit_coeffs(11).*exp(-fit_coeffs(12).*x);
elseif strcmp(fit_type,'polynomial')
    fit_coeffs = B_LIST_polynomial{1,1}{1,12};
    r_fit = fit_coeffs(12).*x.^0 + fit_coeffs(11).*x.^1 + fit_coeffs(10).*x.^2 + fit_coeffs(9).*x.^3 + fit_coeffs(8).*x.^4 + fit_coeffs(7).*x.^5 + ...
        fit_coeffs(6).*x.^6 + fit_coeffs(5).*x.^7 + fit_coeffs(4).*x.^8 + fit_coeffs(3).*x.^9 + fit_coeffs(2).*x.^10 + fit_coeffs(1).*x.^11;
    % fit_coeffs = B_LIST_polynomial{1,1}{1,6};
    % r_fit = fit_coeffs(6).*x.^0 + fit_coeffs(5).*x.^1 + fit_coeffs(4).*x.^2 + fit_coeffs(3).*x.^3 + fit_coeffs(2).*x.^4 + fit_coeffs(1).*x.^5;
elseif strcmp(fit_type,'logistic')
    fit_coeffs = B1(:,1);
    r_fit = fit_coeffs(5) + (fit_coeffs(6)-fit_coeffs(5))./(fit_coeffs(1)+fit_coeffs(3)*exp(-fit_coeffs(4).*x)).^(1/fit_coeffs(2));
end
if TotalRuns == 1
      spiralDeltaV_dt_r_omega0(dt_out,r_fit,omega_list,inertia_ratio_list,1);  
      DeltaV_Optimal = spiralDeltaV_dt_r_omega0(dt_out,yx,omega_list,inertia_ratio_list,0)
      
      for i = 1:num_terms
          switch i
              case 1
                  fit_coeffs = B1(:,end);
                  r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x);
              case 2
                  fit_coeffs = B2(:,end);
                  r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x);
              case 3
                  fit_coeffs = B3(:,end);
                  r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x) + fit_coeffs(5).*exp(-fit_coeffs(6).*x);
              case 4
                  fit_coeffs = B4(:,end);
                  r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x) + fit_coeffs(5).*exp(-fit_coeffs(6).*x) +...
                      fit_coeffs(7).*exp(-fit_coeffs(8).*x);
              case 5
                  fit_coeffs = B5(:,end);
                  r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x) + fit_coeffs(5).*exp(-fit_coeffs(6).*x) +...
                      fit_coeffs(7).*exp(-fit_coeffs(8).*x) + fit_coeffs(9).*exp(-fit_coeffs(10).*x);
              case 6
                  fit_coeffs = B6(:,end);
                  r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x) + fit_coeffs(5).*exp(-fit_coeffs(6).*x) +...
                      fit_coeffs(7).*exp(-fit_coeffs(8).*x) + fit_coeffs(9).*exp(-fit_coeffs(10).*x) + fit_coeffs(11).*exp(-fit_coeffs(12).*x);
          end
          R_FIT(:,i) = fcnfit(i,:);
          DeltaV_Fit(i) = spiralDeltaV_dt_r_omega0(dt_out,fcnfit(i,:),omega_list,inertia_ratio_list,0);
          
      end
      DeltaV_Fit_Normalized = DeltaV_Fit/DeltaV_Optimal;
      if num_terms_best_exponential == 2
          fit_coeffs = [8.37882458523862;0.0551170487169159;1.62117541476138;0.00268646208704073];
          r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x);
          DeltaV_TwoTermBestFit = spiralDeltaV_dt_r_omega0(dt_out,r_fit,omega_list,inertia_ratio_list,0)
      elseif num_terms_best_exponential == 3
          fit_coeffs = [9.04079194195414;0.0509358869930513;2.16063608354774;-0.00805623257712262;-1.20142802550188;-0.0106790005587400];
          r_fit = fit_coeffs(1).*exp(-fit_coeffs(2).*x) + fit_coeffs(3).*exp(-fit_coeffs(4).*x) + fit_coeffs(5).*exp(-fit_coeffs(6).*x);
          DeltaV_ThreeTermBestFit = spiralDeltaV_dt_r_omega0(dt_out,r_fit,omega_list,inertia_ratio_list,0)
      end

      figure()
      plot(x,R_FIT);
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
      
end

% %% Fit an Archimedes Spiral
% if plot_single == 1
%     
%     y = @(b,x) b(1).*(deg2rad(5).*x).^(1./b(2));
%     x = t_out;
%     yx = r;
%     OLS = @(b) sum((y(b,x) - yx).^2);
%     opts = optimset('MaxFunEvals',1E100, 'MaxIter',1E100);
%     B_Archimedes = fminsearch(OLS, [3 2], opts);
%     fcnfitArchimedes = y(B_Archimedes,x);
%     error_in_this_fit = norm(abs(fcnfitArchimedes-yx));
%     fit_error_Archimedes = error_in_this_fit;
%     
%     figure()
%     plot(x, yx, '*-b')
%     hold on
%     plot(x, fcnfitArchimedes, '-r')
%     hold off
%     grid
%     title({'Radial Distance Archimedes Spiral Fit with';['Norm(Error)=',num2str(fit_error_Archimedes)]})
%     xlabel('Time [s]')
%     ylabel('Radial Distance [m]')
%     
%     omega = deg2rad(5);
%     r0 = 10;
%     rf = 1;
%     tf = (1/omega) * ((r0-rf)./B_Archimedes(1)).^B_Archimedes(2);
%     t_archimedes = 0:0.01:tf;
%     r_archimedes = 10-B_Archimedes(1)*(deg2rad(5).*t_archimedes).^(1./B_Archimedes(2));
%     theta_archimedes = omega.*t_archimedes;
%     theta_optimum = omega.*x;
%     figure()
%     polar(theta_optimum,yx,'b'); hold on;
%     polar(theta_archimedes,r_archimedes,'r')
%     legend('Optimum','Archimedes')
%     title('Best-Fit Archimedes Spiral to Optimal Approach')
%     
%     DV_archimedes_fit = spiralDeltaV_dt_r(dt_out,fcnfitArchimedes);
%     
%     percent_difference = 100 * (DV_archimedes_fit/DV_out)
%     
% end
% 
% %%
% if plot_single == 0
%     figure();
%     for i = 1:size(optimal_radial_distribution,1)
%         plot(optimal_time_history(i,:),optimal_radial_distribution(i,:)); hold on;
%     end
%     title('Optimal Radial Distance Profile')
%     xlabel('Time [s]')
%     ylabel('Radius [m]')
%     grid on
% end
% 
% %% Exponential Fit Analysis
% % if strcmp(fit_type,'exponential') || strcmp(fit_type,'exponential_with_r0')
%     B_LIST = B_LIST_exponential;
%     R_LIST = R_LIST_exponential;
%     T_LIST = T_LIST_exponential;
%     DV_LIST = DV_LIST_exponential;
%     OMEGA_LIST = OMEGA_LIST_exponential;
%     RESULTS = RESULTS_exponential;
%     Exponential_Fit_Analysis
% % end
% close all;
% %% Polynomial Fit Analysis
% % if strcmp(fit_type,'polynomial')
%     B_LIST = B_LIST_polynomial;
%     R_LIST = R_LIST_polynomial;
%     T_LIST = T_LIST_polynomial;
%     DV_LIST = DV_LIST_polynomial;
%     OMEGA_LIST = OMEGA_LIST_polynomial;
%     RESULTS = RESULTS_polynomial;
%     Polynomial_Fit_Analysis
% % end
% 
% 
% %% Plot Both Exponential and Polynomial Fits Together
% if size(RESULTS_polynomial,1) == 1
%     
%     % Plot based on number of terms
%    for i = 1:size(RESULTS_polynomial,2) %number of terms
%       figure()
%       plot(optimal_time_history,optimal_radial_distribution,'b-*'); hold on;
%       c_exp = RESULTS_exponential{i}(8:end-3);
%       c_poly = flipud(RESULTS_polynomial{i}(8:end-3));
%       exp_fit_radius = zeros(1,length(optimal_time_history));
%       poly_fit_radius = zeros(1,length(optimal_time_history));
%       for j = 1:2:length(c_exp)
%            exp_fit_radius = exp_fit_radius + c_exp(j)*exp(optimal_time_history*-c_exp(j+1));
%       end
%       plot(optimal_time_history,exp_fit_radius,'r-');
%       for j = 1:length(c_poly)
%            poly_fit_radius = poly_fit_radius + c_poly(j)*optimal_time_history.^(j-1);
%       end
%       plot(optimal_time_history,poly_fit_radius,'k:');
%       legend('Optimal Trajectory','Exponential Fit','Polynomial Fit')
%       xlabel('Time [s]')
%       ylabel('Radial Distance [m]')
%       grid on;
%       title({['Exponential and Polynomial Fits to'];['Optimal Trajectory with ',num2str(i),' Terms']});
%    end
%    
%    
%    % Plot based on number of coefficients
%    for i = 1:size(RESULTS_polynomial,2)/2 %number of coefficients
%       figure()
%       plot(optimal_time_history,optimal_radial_distribution,'b-*'); hold on;
%       c_exp = RESULTS_exponential{i}(8:end-3);
%       c_poly = flipud(RESULTS_polynomial{2*i}(8:end-3));
%       exp_fit_radius = zeros(1,length(optimal_time_history));
%       poly_fit_radius = zeros(1,length(optimal_time_history));
%       for j = 1:2:length(c_exp)
%            exp_fit_radius = exp_fit_radius + c_exp(j)*exp(optimal_time_history*-c_exp(j+1));
%       end
%       plot(optimal_time_history,exp_fit_radius,'r-');
%       for j = 1:length(c_poly)
%            poly_fit_radius = poly_fit_radius + c_poly(j)*optimal_time_history.^(j-1);
%       end
%       plot(optimal_time_history,poly_fit_radius,'k:');
%       legend('Optimal Trajectory','Exponential Fit','Polynomial Fit')
%       xlabel('Time [s]')
%       ylabel('Radial Distance [m]')
%       grid on;
%       title({['Exponential and Polynomial Fits to'];['Optimal Trajectory with ',num2str(2*i),' Coefficients']});
%    end
%    
%    
%     % Plot the fit error over number of coefficients
%     figure()
%     subplot(2,1,1)
%     plot(fit_error_polynomial(:,1),fit_error_polynomial(:,2),'b-'); hold on;
%     plot(2*fit_error_exponential(:,1),fit_error_exponential(:,2),'r-.'); hold on;
%     title('Fit Error Over Number of Coefficients')
%     xlabel('Number of Coefficients') 
%     ylabel('Fit Error')
%     legend('Polynomial Fit','Exponential Fit');
%     subplot(2,1,2)
%     plot(fit_error_polynomial(:,1),fit_error_polynomial(:,2),'b-'); hold on;
%     plot(fit_error_exponential(:,1),fit_error_exponential(:,2),'r-.'); hold on;
%     title('Fit Error Over Number of Terms')
%     xlabel('Number of Terms') 
%     ylabel('Fit Error')
%     legend('Polynomial Fit','Exponential Fit');
% end
% 
