clear all
close all


%% define parameters

dim = 2;

epsilon = 1e-3;
eps_dmaps = 10;
dt = 1e-3;

nsteps = 3000;
nsteps_burst = 100;

%% get initial data

rng(123);
[data_init, ~, data1, ~, t] = simulate_quadratic_linear(nsteps, dt, 1, dt, epsilon);


make_fig;
scatter(data1(:,1),data1(:,2),50,t*dt,'.')
xlabel('Y(1)')
ylabel('Y(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal
print(gcf, '-depsc', 'dat_varyburst');

%%

fig_idx = 1;

for dt_burst=[1e-8 1e-4 1e-2]
    
    rng(123);
    [~, data_burst_init, data1, data1_burst, ~] = simulate_quadratic_linear(nsteps, dt, nsteps_burst, dt_burst, epsilon, data_init);
    
    make_fig;
    %     scatter(data1(:,1),data1(:,2),50,t,'.')
    plot(data1(:,1),data1(:,2),'.b')
    hold on
    plot(data1_burst(:,1,1000),data1_burst(:,2,1000),'.r')
    xlabel('Y(1)')
    ylabel('Y(2)')
    %     h = colorbar('peer',gca);
    %     set(get(h,'xlabel'),'String', 't');
    axis equal
    print(gcf, '-depsc', sprintf('data_withburst_%d', fig_idx));
    
    
    % DMAPS
    
    W = squareform(pdist(data_init)).^2;
    [V_init_dmaps, D_init_dmaps] = dmaps(W, eps_dmaps, 10);
    
    
    % compute covariances
    [inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);
    inv_c1 = inv_c1 * dt_burst;
    
    % NIV
    eps_niv = eps_dmaps;
    [V1_niv, D1_niv, ~, Dis_data1] = NIV_return_dist(data1, inv_c1, eps_niv, 50, 0);
    
    [~, fast_idx] = max(corr(V1_niv, data_init(:,2)));
    
    make_fig;
    plot(data_init(:,2), V1_niv(:,fast_idx),'.')
    ylabel(sprintf('NIV from Y data (\\psi_{%d})', fast_idx))
    xlabel('fast variable')
    print(gcf, '-depsc', sprintf('fast_var_corr_%d', fig_idx));
    
    
    fast_idx
    
    fig_idx = fig_idx + 1;
    
end

