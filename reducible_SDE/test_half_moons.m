% this code simulates the "half-moons" from "Detecting intrinsic slow variables in stochastic
% dynamical systems by anisotropic diffusion maps", Singer, Erban,
% Kevrekidis, and Coifman, PNAS

clear all
close all

%% define parameters
a1 = 1e-3;
a2 = 1e-3;
a3 = 1e-1;
a4 = 1e-1;

% initial conditions
data0 = [0; 1];

% time to integrate
<<<<<<< Updated upstream
tmax = 3e3;
=======
tmax = 2e3;
>>>>>>> Stashed changes

%% do simulations

drift = @(t, x) [a1; a3*(1-x(2))];
diffn = @(t, x) [a2 0; 0 a4];
    
dt = 1;
nsteps_per_step = 1000;

SDE = sde(drift, diffn, 'StartState', data0);
nPeriods = ceil(tmax / dt);
[data, t, Z] = SDE.simulate(nPeriods, 'DeltaTime', dt, 'NSTEPS', nsteps_per_step);

u = data(:,1);
v = data(:,2);

%f = @(x) [x(:,2).*cos(x(:,1)+x(:,2)-1) x(:,2).*sin(x(:,1)+x(:,2)-1)];
f = @(x) x;

fdata = f(data);
x = fdata(:,1);
y = fdata(:,2);

%% plot results

% plot variables, color by time
figure;
scatter(x, y ,50,t,'.')
xlabel('x')
ylabel('y')
title('colored by time')

%% calculate NIV

neigs = 60;
dim = 2;

dt2_range = [1e-4 1e-2 1];

figure;
plot_idx = 1;
m_subfig = length(dt2_range);
n_subfig = 3;

fast_coord = [];
lambda_ratio = [];

for dt2=dt2_range
    knn_NIV = 50;
    nsteps_per_step2 = max([1 round(dt2/1e-3)]);

    data_snippets = zeros(2*knn_NIV+1, size(data, 2), size(data, 1));
    
    if dt2 < dt
        for i=1:size(data, 1) 
            SDE = sde(drift, diffn, 'StartState', data(i,:)');
            [tmp_data, ~, ~] = SDE.simulate(2*knn_NIV, 'DeltaTime', dt2, 'NSTEPS', nsteps_per_step2);
            data_snippets(:,:,i) = f(tmp_data);
        end
    else
        stride = round(dt2/dt);
        for i=1:size(data,1)
            i1 = i;
            i2 = i+2*knn_NIV*stride;
            if i2 > size(data,1)
                shift = i2 - size(data,1);
                i2 = i2 - shift;
                i1 = i1 - shift;
            end
            data_snippets(:,:,i) = fdata(i1:stride:i2,:);
        end
    end
    neigh_idx = floor(size(data,1)/2);

    subplot(m_subfig, n_subfig,plot_idx)
    plot_idx = plot_idx + 1;
    %figure;
    scatter(fdata(:,1),fdata(:,2),50,t,'.')
    hold on
    plot(data_snippets(:, 1, neigh_idx), data_snippets(:, 2, neigh_idx), '.k')
    xlabel('x')
    ylabel('y')
    title('colored by time')
    
    [inv_c, subsampled_data, ranks] = covariances2(data_snippets, dim);
    
    tmp_lambda_ratio = zeros(size(data,1), 1);
    for i=1:size(data,1)
        [u, s, v] = svd(inv_c(:,:,i));
        tmp_lambda_ratio(i) = s(1,1)/s(2,2);
    end
    lambda_ratio = [lambda_ratio tmp_lambda_ratio];
    
    if size(subsampled_data,1) > 4000
        disp('Too much data; NIV will fail')
        return
    end
    [V, D] = NIV(subsampled_data, inv_c, 0, neigs, 0);
    
    % plot NIV
    % define indices for uncorrelated NIV components
    psi_idx = zeros(2,1);
    [~, psi_idx(1)] = max(abs(corr(V, data(:, 1))));
    [~, psi_idx(2)] = max(abs(corr(V, data(:, 2))));
    
    fast_coord = [fast_coord psi_idx(2)];
    
    % show correlation of first NIV with slow variable
    %figure;
    subplot(m_subfig, n_subfig,plot_idx)
    plot_idx = plot_idx + 1;
    plot(V(:,psi_idx(1)),data(:,1),'.')
    xlabel(sprintf('\\psi_{%d}',psi_idx(1)))
    ylabel('u')
    
    % show correlation of second NIV with fast variable
    %figure;
    subplot(m_subfig, n_subfig,plot_idx)
    plot_idx = plot_idx + 1;
    plot(V(:,psi_idx(2)),data(:,2),'.')
    xlabel(sprintf('\\psi_{%d}',psi_idx(2)))
    ylabel('v')
    
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
