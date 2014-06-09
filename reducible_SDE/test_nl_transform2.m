clear all
close all


%%

dim = 2;

a = 1;
epsilon = 0.01;

DriftFn = @(t, x) [a; -x(2)^15/epsilon];
DiffnFn = @(t, x) [0.1 0; 0 1/sqrt(epsilon)];

dt = 0.5e-3;
dt_step = 1e-4;
nsteps = 4000;

nsteps_per_step = dt / dt_step;
if mod(nsteps_per_step, 1) ~= 0
    disp('Step dt is not a multiple of dt');
    return;
end

data_start = [0; 0];

%%

rng(123);

SDE = sde(DriftFn, DiffnFn, 'StartState', data_start);
[data_init, t, Z] = SDE.simulate(nsteps, 'DeltaTime', dt, 'NSTEPS', nsteps_per_step);
data_init = data_init(2:end, :);
t = t(2:end);

%%

figure;
scatter(data_init(:,1),data_init(:,2),50,t)
axis equal 

%%

nsteps_burst = 200;
dt_burst = dt/100;
dt_burst_step = dt_burst;

nsteps_burst_per_step = dt_burst / dt_burst_step;
if mod(nsteps_burst_per_step, 1) ~= 0
    disp('Step dt_burst is not a multiple of dt_burst');
    return;
end

data1 = zeros(nsteps_burst, dim, nsteps);
for i=1:nsteps
    SDE = sde(DriftFn, DiffnFn, 'StartState', data_init(i, :)');
    [data_tmp, ~, Z] = SDE.simulate(1, 'DeltaTime', dt_burst, 'NSTEPS', nsteps_burst_per_step, 'ntrials', nsteps_burst);
    data1(:, :, i) = squeeze(data_tmp(end, :, :))';
end

figure;
scatter(data1(1,1,:),data1(1,2,:),50,t)
axis equal 

figure;
scatter(data1(1,1,:),data1(1,2,:))
hold on
plot(data1(:, 1, 1000), data1(:,2,1000),'.r')


%%


% f = @(x) [(x(:,2)+5).*cos(x(:,1)+transform_curvature*x(:,2)) (x(:,2)+5).*sin(x(:,1)+transform_curvature*x(:,2))];

% R = max(abs(data1(:,2)))+0.01;
%f = @(x) [x(:,1)+sqrt(R^2-(x(:,2)).^2) x(:,2)];

freq = 0.5;
scale = 1;
f = @(x) [x(:,1)+scale*sin(x(:,2)/freq) x(:,2)];

data2 = zeros(size(data1));
for i=1:nsteps
    data2(:, :, i) = f(data1(:, :, i));
end

figure;
scatter(data2(1,1,:),data2(1,2,:),50,t)
axis equal 

%% DMAPS

W = squareform(pdist(squeeze(data1(1,:,:))')).^2;
eps = median(W(:));

[V1_dmaps, D1_dmaps] = dmaps(W, eps, 10);

figure;
scatter(data1(1,1,:),data1(1,2,:),50,V1_dmaps(:,2))
axis equal


%%
W = squareform(pdist(squeeze(data2(1,:,:))')).^2;
eps = median(W(:));

[V2_dmaps, D2_dmaps] = dmaps(W, eps, 10);

figure;
scatter(data2(1,1,:),data2(1,2,:),50,V2_dmaps(:,2))
axis equal

figure;
scatter(data1(1,1,:),data1(1,2,:),50,V2_dmaps(:,2))
axis equal


figure;
plot(V1_dmaps(:,2), V2_dmaps(:,2),'.')


corr(V1_dmaps(:,2), V2_dmaps(:,2))


%%

[inv_c1, new_data1, ranks1] = covariances2(data1, dim);
[inv_c2, new_data2, ranks2] = covariances2(data2, dim);

%%

[V1_niv, D1_niv, eps] = NIV(new_data1, inv_c1, 0, 10, 0);


figure;
scatter(new_data1(:,1),new_data1(:,2),50,V1_niv(:,2))
axis equal

%%

[V2_niv, D2_niv, eps] = NIV(new_data2, inv_c2, 1e6, 10, 0);

figure;
scatter(new_data2(:,1),new_data2(:,2),50,V2_niv(:,2))
axis equal


figure;
scatter(new_data1(:,1),new_data1(:,2),50,V2_niv(:,2))
axis equal

%%


figure;
plot(V1_niv(:,2), V2_niv(:,2),'.')


corr(V1_niv(:,2), V2_niv(:,2))
