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
tmax = 10;
% Gillespie steps to take in each (saved) step of the integration
steps_per_it = 1000;

%% do Gillespie simulations
[t, n] = Gillespie(n0, tmax, steps_per_it);

%% plot results

% plot fast and slow variables after concentrations have relaxed
time_idx = find(t > 0.01, 1, 'first');
figure;
plot(n(time_idx:end,slow_idx),n(time_idx:end,fast_idx),'.')
xlabel(slow_idx_label)
ylabel(fast_idx_label)

% plot variables, color by time
figure;
scatter(n(time_idx:end,slow_idx),n(time_idx:end,fast_idx),50,t(time_idx:end),'.')
xlabel(slow_idx_label)
ylabel(fast_idx_label)
title('colored by time')

%% calculate NIV
stride = 3;
neigs = 20;

inv_c = covariances(n(time_idx:stride:end,:), 50, 2);
[V, D] = NLICA(n(time_idx:stride:end,:), inv_c, neigs);

%% plot NIV

% define indices for uncorrelated NIV components
psi_idx = [2, 10];

% show 2D NIV embedding
figure; 
plot(V(:,psi_idx(1)),V(:,psi_idx(2)),'.')
xlabel(sprintf('\\psi_{%d}',psi_idx(1)))
ylabel(sprintf('\\psi_{%d}',psi_idx(2)))

% show 2D NIV embedding colored by time
figure; 
scatter(V(:,psi_idx(1)),V(:,psi_idx(2)),50,t(time_idx:stride:end),'.')
xlabel(sprintf('\\psi_{%d}',psi_idx(1)))
ylabel(sprintf('\\psi_{%d}',psi_idx(2)))
title('colored by time')

% show correlation of first NIV with slow variable
figure; 
plot(V(:,psi_idx(1)),n(time_idx:stride:end,slow_idx),'.')
xlabel(sprintf('\\psi_{%d}',psi_idx(1)))
ylabel(slow_idx_label)

% show correlation of second NIV with fast variable
figure; 
plot(V(:,psi_idx(2)),n(time_idx:stride:end,fast_idx),'.')
xlabel(sprintf('\\psi_{%d}',psi_idx(2)))
ylabel(fast_idx_label)
