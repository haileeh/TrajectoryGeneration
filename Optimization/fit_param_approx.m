function [fcnfit,fit_error,B_list]=fit_param_approx(optimal_radial_distribution,x,...
    r0, rf, fit_type, NUM_TERMS,num_dt)

fcnfit = zeros(NUM_TERMS,num_dt+1);

num_evals_iters = [1E10, 1E100];
n = length(num_evals_iters);
B1 = zeros(n,2); B2 = zeros(2*n,2); B3 = zeros(3*n,2); B4 = zeros(4*n,2);...
    B5 = zeros(5*n,2); B6 = zeros(6*n,2);
for num_terms = 1:NUM_TERMS
    % For exponential fit
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
    
    % Objective function for fmincon in block below
    OLS = @(b) sum((y(b,x) - optimal_radial_distribution').^2);
    % OLS is finding the least squares wrt optimized trajectory (r_fit - r_opt)
    
    for i = 1:n %num_evals_iters = [1E10, 1E100]
        opts = optimoptions(@fmincon,'ConstraintTolerance',1e-9,...
            'MaxFunctionEvaluations',num_evals_iters(i),'MaxIterations',num_evals_iters(i),...
            'FunctionTolerance',1e-7,'Display','Off');
        switch num_terms
            case 1
                nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                B = fmincon(OLS, [r0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(2,1),[],nonlcon,opts);
                % fmincon(objective, x0, A, b, Aeq, beq, lower bounds,
                % upper bounds, nonlcon, opts)
                B1(:,i) = B; %[B1,B];
            case 2
                nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                B = fmincon(OLS, [B_list{1};0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(4,1),[],nonlcon,opts);
                B2(:,i) = B;%[B2,B];
            case 3
                nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                B = fmincon(OLS, [B_list{2};0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(6,1),[],nonlcon,opts);
                B3(:,i) = B; %[B3,B];
            case 4
                nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                B = fmincon(OLS, [B_list{3};0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(8,1),[],nonlcon,opts);
                B4(:,i) = B;%[B4,B];
            case 5
                nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                B = fmincon(OLS, [B_list{4};0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(10,1),[],nonlcon,opts);
                B5(:,i) = B;%[B5,B];
            case 6
                nonlcon = @(b)TrajectoryOptimizerVariableStudy_nonlconfunc_B(x,b,r0,rf,fit_type);
                B = fmincon(OLS, [B_list{5};0;0],[], [],repmat([1,0],1,num_terms),r0,-10*ones(12,1),[],nonlcon,opts);
                B6(:,i) = B;%[B6,B];
        end
        
        % Error checking
        fcnfit(num_terms,:) = y(B,x); % Calculate function with estimated parameters
        error_in_this_fit = norm(abs(fcnfit(num_terms,:)-optimal_radial_distribution));
        fit_error(num_terms,:) = [num_terms, error_in_this_fit, num_evals_iters(i)];
        
        %RESULTS_exponential{run_num_exponential,num_terms} = [r0;rf;omega;Tend;error_in_this_fit;B;inertia_ratio];
        B_list{num_terms} = B;
        %fit_error_exponential = [fit_error_exponential;fit_error(num_terms,:)];
        
    end %loop over num_evals_iters
    
    figure()
    plot(x, optimal_radial_distribution, '*-b'); hold on
    plot(x, fcnfit(num_terms,:), '-r'); hold off
    grid
    title({'Radial Distance Exponential Fit with';[num2str(num_terms),' Terms; Norm(Error)=',num2str(norm(fcnfit(num_terms,:)-optimal_radial_distribution))]})
    xlabel('Time [s]')
    ylabel('Radial Distance [m]')
    legend('Optimal','Exponential')
    
end % loop over num_terms

figure();
plot(fit_error(:,1),fit_error(:,2),'-*k')
title('Fit Error over Number of Exponential Terms')
xlabel('Number of Terms')
ylabel('Norm of Error')