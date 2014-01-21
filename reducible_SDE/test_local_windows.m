% this code simulates the detailed example from Section 5.3 of
% Contou-Carrere et al., "Model reduction of multi-scale chemical Langevin
% equations," Systems & Control Letters, 60, 2011.

clear all
close all

%% define parameters
% rate constants, etc. defined in this file
parameters;

% initial conditions
n0 = [3000;
    2500;
    1500;
    2000;
    3000];

% slow and fast variable indices and labels
slow_idx = 4;
fast_idx = 1;
slow_idx_label = 'X_D';
fast_idx_label = 'X_A';

% time to integrate
tmax = 2;
%tmax = 0.5;

%% do simulations

% Gillespie steps to take in each (saved) step of the integration
%steps_per_it = 1000;
%[t, n] = Gillespie(n0, tmax, steps_per_it);

dt = 1e-5;
nsteps_per_step = 1;
[t, n] = SDE_simulation(n0, tmax, dt, nsteps_per_step);

%% plot results

% % plot fast and slow variables after concentrations have relaxed
time_idx = find(t > 0.01, 1, 'first');
n_all = n;
t_all = t;
n = n_all(time_idx:end,:);
t = t_all(time_idx:end,:);

% figure;
% plot(n(:,slow_idx),n(:,fast_idx),'.')
% xlabel(slow_idx_label)
% ylabel(fast_idx_label)
% 
% % plot variables, color by time
% figure;
% scatter(n(:,slow_idx),n(:,fast_idx),50,t,'.')
% xlabel(slow_idx_label)
% ylabel(fast_idx_label)
% title('colored by time')

%% calculate NIV
% figure;
% plot_idx = 1;
% m_subfig = 3;
% n_subfig = 3;

stride = 50;

%for knn_NIV=[5 50 200]*stride
for knn_NIV=20
    neigs = 60;
    dim = 2;
    
    neigh_idx = floor(size(n,1)/2);

    %subplot(m_subfig, n_subfig,plot_idx)
    %plot_idx = plot_idx + 1;
    figure;
    scatter(n(:,slow_idx),n(:,fast_idx),50,t,'.')
    hold on
    plot(n(neigh_idx-knn_NIV:neigh_idx+knn_NIV, slow_idx), n(neigh_idx-knn_NIV:neigh_idx+knn_NIV, fast_idx), '.k')
    xlabel(slow_idx_label)
    ylabel(fast_idx_label)
    title('colored by time')
    
    
    [inv_c, subsampled_n, ranks] = covariances(n, knn_NIV, dim, stride);
    
%     subplot(m_subfig, n_subfig,plot_idx)
%     plot_idx = plot_idx + 1;
%     plot(ranks, '.')
    
    if size(subsampled_n,1) > 4000
        disp('Too much data; NIV will fail')
        return
    end
    [V, D] = NIV(subsampled_n, inv_c, 0, neigs, 0);
    
    % plot NIV
    % define indices for uncorrelated NIV components
    psi_idx = zeros(2,1);
    [~, psi_idx(1)] = max(abs(corr(V, subsampled_n(:, slow_idx))));
    [~, psi_idx(2)] = max(abs(corr(V, subsampled_n(:, fast_idx))));
    %psi_idx = [2, 6];
    
    % show correlation of first NIV with slow variable
    figure;
    %subplot(m_subfig, n_subfig,plot_idx)
    %plot_idx = plot_idx + 1;
    plot(V(:,psi_idx(1)),subsampled_n(:,slow_idx),'.')
    xlabel(sprintf('\\psi_{%d}',psi_idx(1)))
    ylabel(slow_idx_label)
    
    % show correlation of second NIV with fast variable
    figure;
    %subplot(m_subfig, n_subfig,plot_idx)
    %plot_idx = plot_idx + 1;
    plot(V(:,psi_idx(2)),subsampled_n(:,fast_idx),'.')
    xlabel(sprintf('\\psi_{%d}',psi_idx(2)))
    ylabel(fast_idx_label)
    
end
