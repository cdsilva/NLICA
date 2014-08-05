clear all
close all


%% define parameters

epsilon = 1e-3;
eps_dmaps_all = [0.01 10 0.01];
dt_burst_all = [1e-9 1e-9 1e-5];

dim = 2;

dt = 1e-4;

nsteps = 3000;
nsteps_burst = 200;

dy_axis_lim = [1e-2 1e1 1e-10 1e5];


for i = 1:length(eps_dmaps_all)
    
    fig_idx = i;
    eps_dmaps = eps_dmaps_all(i);
    dt_burst = dt_burst_all(i);
    
    make_fig;
    deltas = logspace(-3, 1, 100);
    loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2)
    hold on
    loglog(deltas, (38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2, '-g')
    loglog(deltas, 10*deltas.^4, '-r')
    loglog(sqrt(eps_dmaps/(0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2)))), eps_dmaps, 'ok')
    % legend('linear approximation distance','error from covariance estimation','error from Taylor expansion','location','best')
    xlabel('$\| Y_2 - Y_1 \|$','interpreter','latex')
    ylabel('$\| X_2 - X_1 \|^2_{est}$', 'interpreter','latex')
    axis(dy_axis_lim)
    print(gcf, '-depsc', sprintf('changing_parameters_error_terms_%d', fig_idx));
    
    
    %% simulate SDE
    
    rng(123);
    [data_init, data_burst_init, data1, data1_burst, t] = simulate_quadratic(nsteps, dt, nsteps_burst, dt_burst, epsilon);
    
    
    % make_fig;
    % scatter(data_init(:,1),data_init(:,2),50,t,'.')
    % xlabel('X(1)')
    % ylabel('X(2)')
    % h = colorbar('peer',gca);
    % set(get(h,'xlabel'),'String', 't');
    % axis equal
    % % print(gcf, '-depsc', sprintf('orig_data_%d', fig_idx));
    %
    % make_fig;
    % scatter(data1(:,1),data1(:,2),50,t,'.')
    % xlabel('Y(1)')
    % ylabel('Y(2)')
    % h = colorbar('peer',gca);
    % set(get(h,'xlabel'),'String', 't');
    % axis equal
    % % print(gcf, '-depsc', sprintf('function_data_%d', fig_idx));
    
    
    %%
%     [inv_c_init, new_data_init, ranks_init] = covariances2(data_burst_init, dim);
%     inv_c_init = inv_c_init * dt_burst;
    
    [inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);    
    inv_c1 = inv_c1 * dt_burst;
    
    %%
    eps_niv = eps_dmaps;
    [V1_niv, D1_niv, ~, Dis_data1] = NIV_return_dist(data1, inv_c1, eps_niv, 10, 0);
    
    make_fig;
    plot(data_init(:,1), V1_niv(:,2), '.')
    axis([-inf inf -inf inf])
    xlabel('X(1)')
    ylabel('\phi_1')
    print(gcf, '-depsc', sprintf('changing_parameters_NIV_corr_%d', fig_idx));
    
end
