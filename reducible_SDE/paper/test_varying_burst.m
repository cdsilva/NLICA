clear all
close all


%% define parameters

dim = 2;

epsilon = 1e-3;


%%

for dt_burst=[1e-9 1e-4 1e-2]
    
    eps_dmaps = 1;
    % dt_burst = 1e-9;
    dt = 1e-3;
    
    
    nsteps = 3000;
    nsteps_burst = 100;
    
    
    rng(123);
    [data_init, data_burst_init, data1, data1_burst, t] = simulate_quadratic_linear(nsteps, dt, nsteps_burst, dt_burst, epsilon);
    
%     
%     make_fig;
%     scatter(data_init(:,1),data_init(:,2),50,t,'.')
%     xlabel('X(1)')
%     ylabel('X(2)')
%     h = colorbar('peer',gca);
%     set(get(h,'xlabel'),'String', 't');
%     axis equal
    
%     make_fig;
%     scatter(data1(:,1),data1(:,2),50,t,'.')
%     xlabel('Y(1)')
%     ylabel('Y(2)')
%     h = colorbar('peer',gca);
%     set(get(h,'xlabel'),'String', 't');
%     axis equal
%     
    make_fig;
%     scatter(data1(:,1),data1(:,2),50,t,'.')
    plot(data1(:,1),data1(:,2),'.r')
    hold on
    plot(data1_burst(:,1,1000),data1_burst(:,2,1000),'.k')
    xlabel('Y(1)')
    ylabel('Y(2)')
    h = colorbar('peer',gca);
    set(get(h,'xlabel'),'String', 't');
    axis equal
    
    %% DMAPS
    
    W = squareform(pdist(data_init)).^2;
    [V_init_dmaps, D_init_dmaps] = dmaps(W, eps_dmaps, 10);
    
    
    %%
    
    [inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);
    inv_c1 = inv_c1 * dt_burst;
    
    %%
    
%     eps_niv = eps_dmaps;
    eps_niv = 0;
    
    [V1_niv, D1_niv, ~, Dis_data1] = NIV_return_dist(data1, inv_c1, eps_niv, 50, 0);
    
    
%     make_fig;
%     scatter(data_init(:,1),data_init(:,2),50,V1_niv(:,2),'.')
%     xlabel('x')
%     ylabel('y')
%     h = colorbar('peer',gca);
%     set(get(h,'xlabel'),'String', '\phi_1');
%     axis equal
    
%     make_fig;
%     plot(V1_niv(:,2), V_init_dmaps(:,2),'.')
%     xlabel('NIV from Y data')
%     ylabel('DMAPS from X data')
%     axis equal
% 
%     make_fig;
%     plot(data_init(:,1), V1_niv(:,2),'.')
%     ylabel('NIV from Y data')
%     xlabel('slow variable')

    [~, fast_idx] = max(corr(V1_niv, data_init(:,2)));
    
    make_fig;
    plot(data_init(:,2), V1_niv(:,fast_idx),'.')
    ylabel(sprintf('NIV from Y data \\psi_{%d}', fast_idx))
    xlabel('fast variable')
    
    fast_idx

end

