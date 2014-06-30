clear all
close all


%% define parameters

dim = 2;

a = 3;
epsilon = 0.001;
eps_dmaps = 0.02;

DriftFn = @(t, x) [a; -x(2)/epsilon];
DiffnFn = @(t, x) [1 0; 0 1];
f1 = @(x) [x(:,1)+x(:,2).^2/epsilon x(:,2)/sqrt(epsilon)];

dt = 1e-3;
% dt_burst = 1e-8;
dt_burst = 1e-5;


%%

figure;
deltas = linspace(1e-4, 1, 100);
semilogy(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2)
hold on
semilogy(deltas, 10*deltas.^4, '-r')
semilogy(deltas, (38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2, '-g')
semilogy(sqrt(eps_dmaps/(0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2)))), eps_dmaps, 'ok')
legend('true distance','error from Taylor expansion','error from covariance estimation','location','best')



%%
min_dt = 1e-4;
dt_step = min(dt, min_dt);
nsteps = 3000;

nsteps_per_step = dt / dt_step;
if mod(nsteps_per_step, 1) ~= 0
    disp('Step dt is not a multiple of dt');
    return;
end

data_start = [0; 0];

%% simulate SDE

rng(123);

SDE = sde(DriftFn, DiffnFn, 'StartState', data_start);
[data_init, t, Z] = SDE.simulate(nsteps, 'DeltaTime', dt, 'NSTEPS', nsteps_per_step);
data_init = data_init(2:end, :);
t = t(2:end);

% idx = find(data_init(:, 2).^2 < epsilon/2);
% data_init = data_init(idx, :);
% t = t(idx, :);
% nsteps = length(idx);

%% simulate bursts

nsteps_burst = 200;

dt_burst_step = min(dt_burst, min_dt);

nsteps_burst_per_step = dt_burst / dt_burst_step;
if mod(nsteps_burst_per_step, 1) ~= 0
    disp('Step dt_burst is not a multiple of dt_burst');
    return;
end

data_burst_init = zeros(nsteps_burst, dim, nsteps);
for i=1:nsteps
    SDE = sde(DriftFn, DiffnFn, 'StartState', data_init(i, :)');
    [data_tmp, ~, Z] = SDE.simulate(1, 'DeltaTime', dt_burst, 'NSTEPS', nsteps_burst_per_step, 'ntrials', nsteps_burst);
    data_burst_init(:, :, i) = squeeze(data_tmp(end, :, :))';
end

%% apply functions


data1 = f1(data_init);

data1_burst = zeros(size(data_burst_init));
for i=1:nsteps
    data1_burst(:, :, i) = f1(data_burst_init(:, :, i));
end

%%

make_fig;
scatter(data_init(:,1),data_init(:,2),50,t,'.')
xlabel('X^1')
ylabel('X^2')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal

make_fig;
scatter(data1(:,1),data1(:,2),50,t,'.')
xlabel('Y^1')
ylabel('Y^2')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal




%% DMAPS

W = squareform(pdist(data_init)).^2;
[V_init_dmaps, D_init_dmaps] = dmaps(W, eps_dmaps, 10);

make_fig;
scatter(data_init(:,1),data_init(:,2),50,V_init_dmaps(:,2),'.')
xlabel('x')
ylabel('y')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal


%%
[inv_c_init, new_data_init, ranks_init] = covariances2(data_burst_init, dim);

[inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);

inv_c_init = inv_c_init * dt_burst;
inv_c1 = inv_c1 * dt_burst;

%%

eps_niv = eps_dmaps;
[V_init_niv, D_init_niv, ~, Dis_data_init] = NIV_return_dist(data_init, inv_c_init, eps_niv, 10, 0);


make_fig;
scatter(data_init(:,1),data_init(:,2),50,V_init_niv(:,2),'.')
xlabel('x')
ylabel('y')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal

make_fig;
plot(V_init_niv(:,2), V_init_dmaps(:,2),'.')
xlabel('NIV from fast-slow system')
ylabel('DMAPS from rescaled system')
axis equal

%%
eps_niv = eps_dmaps;
[V1_niv, D1_niv, ~, Dis_data1] = NIV_return_dist(data1, inv_c1, eps_niv, 10, 0);


make_fig;
scatter(data_init(:,1),data_init(:,2),50,V1_niv(:,2),'.')
xlabel('x')
ylabel('y')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal

make_fig;
plot(V1_niv(:,2), V_init_dmaps(:,2),'.')
xlabel('NIV from fast-slow system')
ylabel('DMAPS from rescaled system')
axis equal
