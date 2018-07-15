function num_spy = spy_feasibility_existence(TOTAL_POINTS)
% Function to create a spy plot of when a feasibility plot exists if there
% is a constraint that downselects from the entire set of TOTAL_POINTS as
% output from running the feasibility map code for multiple rf and omega


TOTAL_POINTS(7,:) = floor((10^6*TOTAL_POINTS(7,:)))/10^6; %rounding

omega_list = unique(TOTAL_POINTS(7,:));
num_omega = length(omega_list);
rf_list = unique(TOTAL_POINTS(8,:));
num_rf = length(rf_list);

spy_plot_matrix = zeros(num_omega,num_rf);
for omega_iter = 1:num_omega
    for rf_iter = 1:num_rf
        [~,idx] = find(TOTAL_POINTS(7,:)==omega_list(omega_iter) & (TOTAL_POINTS(8,:)==rf_list(rf_iter)));
        [~,idx_valid] = find(TOTAL_POINTS(3,idx) < 1);
        if idx_valid > 0
           spy_plot_matrix(omega_iter,rf_iter)= 1; 
        end
    end
end

assignin('base','spy_plot_matrix',spy_plot_matrix);

figure(99988)
spy(spy_plot_matrix)
set(gca,'XTick',1:num_rf)
set(gca,'XTickLabel',rf_list)
set(gca,'YTick',1:num_omega)
set(gca,'YTickLabel',round(rad2deg(omega_list),2))
xlabel('Final Radius [m]')
ylabel('Rotation Rate [deg/s]')
set(gca,'FontSize',20)
end