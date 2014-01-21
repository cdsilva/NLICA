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
tmax = 3;

%% do simulations

% Gillespie steps to take in each (saved) step of the integration
%steps_per_it = 1000;
%[t, n] = Gillespie(n0, tmax, steps_per_it);

dt = 1e-3;
nsteps_per_step = 1000;
[t, n] = SDE_simulation(n0, tmax, dt, nsteps_per_step);

%% plot results

% % plot fast and slow variables after concentrations have relaxed
time_idx = find(t > 0.01, 1, 'first');
n_all = n;
t_all = t;
n = n_all(time_idx:end,:);
t = t_all(time_idx:end,:);

figure;
plot(n(:,slow_idx),n(:,fast_idx),'.')
xlabel(slow_idx_label)
ylabel(fast_idx_label)

% plot variables, color by time
figure;
scatter(n(:,slow_idx),n(:,fast_idx),50,t,'.')
xlabel(slow_idx_label)
ylabel(fast_idx_label)
title('colored by time')

%% calculate NIV
figure;
plot_idx = 1;
m_subfig = 2;
n_subfig = 3;

neigs = 60;
dim = 2;

dt2 = 1e-13;
nsteps_per_step2 = 1;

for knn_NIV=[50]
    n_snippets = zeros(2*knn_NIV+1, size(n, 2), size(n, 1));
    
    for i=1:size(n, 1) 
        [~, n_snippets(:,:,i)] = SDE_simulation(n(i,:)', 2*knn_NIV*dt2, dt2, nsteps_per_step2);
    end
    
    neigh_idx = floor(size(n,1)/2);

    subplot(m_subfig, n_subfig,plot_idx)
    plot_idx = plot_idx + 1;
    %figure;
    scatter(n(:,slow_idx),n(:,fast_idx),50,t,'.')
    hold on
    plot(n_snippets(:, slow_idx, neigh_idx), n_snippets(:, fast_idx, neigh_idx), '.k')
    xlabel(slow_idx_label)
    ylabel(fast_idx_label)
    title('colored by time')
    
    
    [inv_c, subsampled_n, ranks] = covariances2(n_snippets, dim);
    
    mean_ratios = zeros(size(n,1),1);
    A_drift = zeros(size(n,1),1);
    D_drift = zeros(size(n,1),1);
    diff_drift = zeros(size(n,1),1);
    for i=1:size(n,1)
        A = DiffusionFn(0, subsampled_n(i,:), k);
        %norm(pinv(A*A') - inv_c(:,:,i)*(dt2*(knn_NIV*2+1)/4)) / norm(pinv(A*A'))
        mean_ratios(i) = mean(mean(pinv(A*A') ./ (inv_c(:,:,i)*(dt2*(knn_NIV*2+1)/(4*3)))));
        %pinv(A*A') ./ inv_c(:,:,i)
        B = DriftFn(0, subsampled_n(i,:), k);
        A_drift(i) = B(1);
        D_drift(i) = B(4);
        tmp_drift1 = DriftFn(0, n_snippets(1,:,i), k);
        tmp_drift2 = DriftFn(0, n_snippets(end,:,i), k);
        diff_drift(i) = norm(tmp_drift1(1)-tmp_drift2(1))/norm(tmp_drift1(1));
    end
    
%     subplot(m_subfig, n_subfig,plot_idx)
%     plot_idx = plot_idx + 1;
%     %figure;
%     scatter(n(:,slow_idx),n(:,fast_idx),50,mean_ratios,'.')
%     xlabel(slow_idx_label)
%     ylabel(fast_idx_label)
%     title('colored by ratio of inverse covariances')
    
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
    %figure;
    subplot(m_subfig, n_subfig,plot_idx)
    plot_idx = plot_idx + 1;
    plot(V(:,psi_idx(1)),subsampled_n(:,slow_idx),'.')
    xlabel(sprintf('\\psi_{%d}',psi_idx(1)))
    ylabel(slow_idx_label)
    
    % show correlation of second NIV with fast variable
    %figure;
    subplot(m_subfig, n_subfig,plot_idx)
    plot_idx = plot_idx + 1;
    plot(V(:,psi_idx(2)),subsampled_n(:,fast_idx),'.')
    xlabel(sprintf('\\psi_{%d}',psi_idx(2)))
    ylabel(fast_idx_label)
    
end
