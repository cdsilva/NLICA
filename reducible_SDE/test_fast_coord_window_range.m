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

%% do simulations

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

neigs = 60;
dim = 2;

%dt2_range = [1e-2 1e-5 1e-8 1e-13];
dt2_range = [1e-3 0.5e-3 1e-4 1e-5 1e-7 1e-10];

figure;
plot_idx = 1;
m_subfig = length(dt2_range);
n_subfig = 3;

fast_coord = [];
lambda_ratio = [];

for dt2=dt2_range
    knn_NIV = 50;
    nsteps_per_step2 = max([1 round(dt2/1e-5)]);

    n_snippets = zeros(2*knn_NIV+1, size(n, 2), size(n, 1));
    
    if dt2 < dt
        for i=1:size(n, 1) 
            [~, n_snippets(:,:,i)] = SDE_simulation(n(i,:)', 2*knn_NIV*dt2, dt2, nsteps_per_step2);
        end
    else
        stride = round(dt2/dt);
        for i=1:size(n,1)
            i1 = i;
            i2 = i+2*knn_NIV*stride;
            if i2 > size(n,1)
                shift = i2 - size(n,1);
                i2 = i2 - shift;
                i1 = i1 - shift;
            end
            n_snippets(:,:,i) = n(i1:stride:i2,:);
        end
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
    
    tmp_lambda_ratio = zeros(size(n,1), 1);
    for i=1:size(n,1)
        [u, s, v] = svd(inv_c(:,:,i));
        tmp_lambda_ratio(i) = s(1,1)/s(2,2);
    end
    lambda_ratio = [lambda_ratio tmp_lambda_ratio];
    
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
    
    fast_coord = [fast_coord psi_idx(2)];
    
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
    
%     subplot(m_subfig, n_subfig,plot_idx)
%     plot_idx = plot_idx + 1;
%     plot(lambda_ratio)
end

figure; 
semilogx(dt2_range,median(lambda_ratio),'.', 'markersize', 12)
xlabel('dt')
ylabel('ratio of \lambda')

figure; 
semilogx(dt2_range,fast_coord,'.', 'markersize', 12)
xlabel('dt')
ylabel('fast index')
