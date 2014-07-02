clear all
close all


%% define parameters
% 
fig_idx = 1;
epsilon = 1e-3;
eps_dmaps = 0.01;
dt_burst = 1e-10;
% 
% fig_idx = 2;
% epsilon = 1e-3;
% eps_dmaps = 10;
% dt_burst = 1e-8;
% 
% fig_idx = 3;
% epsilon = 1e-3;
% eps_dmaps = 0.01;
% dt_burst = 1e-5;
% 
% fig_idx = 4;
% epsilon = 1e-4;
% eps_dmaps = 0.01;
% dt_burst = 1e-10;
% 
% fig_idx = 5;
% epsilon = 1e-4;
% eps_dmaps = 0.01;
% dt_burst = 1e-5;

dim = 2;

a = 3;
DriftFn = @(t, x) [a; -x(2)/epsilon];
DiffnFn = @(t, x) [1 0; 0 1];

f1 = @(x) [x(:,1)+x(:,2).^2/epsilon x(:,2)/sqrt(epsilon)];

dt = 1e-4;
min_dt = 1e-5;

nsteps = 3000;

%%

make_fig;
deltas = logspace(-3, 2, 100);
loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2)
hold on
loglog(deltas, (38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2, '-g')
loglog(deltas, 10*deltas.^4, '-r')
loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2+(38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4, '-k')
% loglog(sqrt(eps_dmaps/(0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2)))), eps_dmaps, 'ok')
% legend('linear approximation distance','error from covariance estimation','error from Taylor expansion','location','best')
xlabel('$\Delta Y$','interpreter','latex')

figure;
loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2+(38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4, '-k')


%%

dt_step = min(dt, min_dt);

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

make_fig;
scatter(data_init(:,1),data_init(:,2),50,t,'.')
xlabel('X(1)')
ylabel('X(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal

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


make_fig;
scatter(data1(:,1),data1(:,2),50,t,'.')
xlabel('Y(1)')
ylabel('Y(2)')
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
xlabel('NIV from Y data')
ylabel('DMAPS from X data')
axis equal

Dis_Y = squareform(pdist(data1));
make_fig;
loglog(Dis_Y(:), Dis_data1(:), '.')
hold on
loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2+(38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4, '-r')

deltas = logspace(-2, 0.5, 100);
curve_means = zeros(size(deltas));
curve_means = curve_means(1:end-1);
[bincounts, ind] = histc(Dis_Y(:), deltas);
for i=1:length(deltas)-1
    curve_means(i) = mean(Dis_data1(ind == i));
end
figure;
loglog(deltas(1:end-1), curve_means, '.b')
hold on
loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2+(38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4, '-r')



%%

make_fig;
deltas = 0.01;
dt_burst = logspace(-10, -2, 100);
loglog(dt_burst, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2*ones(size(dt_burst)))
hold on
loglog(dt_burst, (38 * dt_burst./(epsilon^2 + 2 * dt_burst))*deltas.^2, '-g')
loglog(dt_burst, 10*deltas.^4*ones(size(dt_burst)), '-r')
loglog(dt_burst,  0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2*ones(size(dt_burst))+(38 * dt_burst./(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4*ones(size(dt_burst)), '-k')
xlabel('$\delta t$','interpreter','latex')

%%

nsteps = 10;

dt = 0.0001;
dt_step = min(dt, min_dt);

nsteps_per_step = round(dt / dt_step);
if mod(nsteps_per_step, 1) ~= 0
    disp('Step dt is not a multiple of dt');
    return;
end

data_start = [0; 0];

rng(321);

SDE = sde(DriftFn, DiffnFn, 'StartState', data_start);
[data_tmp, t, Z] = SDE.simulate(1, 'DeltaTime', dt, 'NSTEPS', nsteps_per_step, 'ntrials', nsteps);
data_init = squeeze(data_tmp(end, :, :))';


nsteps_burst = 200;
dt_burst = logspace(-10, -2, 9);

mean_Dis = zeros(size(dt_burst));

for j=1:length(dt_burst)
    dt_burst_step = min(dt_burst(j), min_dt);

    nsteps_burst_per_step = round(dt_burst(j) / dt_burst_step);
%     if mod(nsteps_burst_per_step, 1) ~= 0
%         disp('Step dt_burst is not a multiple of dt_burst');
%         return;
%     end

    data_burst_init = zeros(nsteps_burst, dim, nsteps);
    for i=1:nsteps
        SDE = sde(DriftFn, DiffnFn, 'StartState', data_init(i, :)');
        [data_tmp, ~, Z] = SDE.simulate(1, 'DeltaTime', dt_burst(j), 'NSTEPS', nsteps_burst_per_step, 'ntrials', nsteps_burst);
        data_burst_init(:, :, i) = squeeze(data_tmp(end, :, :))';
    end

    % apply functions
    data1 = f1(data_init);

    data1_burst = zeros(size(data_burst_init));
    for i=1:nsteps
        data1_burst(:, :, i) = f1(data_burst_init(:, :, i));
    end
    

    [inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);
    inv_c1 = inv_c1 * dt_burst(j);

    eps_niv = eps_dmaps;
    [V1_niv, D1_niv, ~, Dis_data1] = NIV_return_dist(data1, inv_c1, eps_niv, 10, 0);

    mean_Dis(j) = mean(Dis_data1(:));
end

deltas = mean(pdist(data1));
figure;
loglog(dt_burst, mean_Dis)
hold on
loglog(dt_burst,  0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2*ones(size(dt_burst))+(38 * dt_burst./(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4*ones(size(dt_burst)), '-k')
