clear all
close all


%% define parameters

fmt = '-djpeg';
res = '-r600';

dim = 2;

a = 3;
epsilon = 0.001;

DriftFn = @(t, x) [a; -x(2)/epsilon];
DiffnFn = @(t, x) [1 0; 0 1/sqrt(epsilon)];

dt = 1e-4;
dt_step = dt/10;
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

%% simulate bursts

nsteps_burst = 100;
dt_burst = dt/10;
dt_burst_step = dt_burst/10;

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

scale = 1*sqrt(epsilon);
f1 = @(x) [x(:,1) scale*x(:,2)];

transform_curvature1 = 2;
transform_curvature2 = 0.5;
f2 = @(x) [(x(:,2)+5).*cos(transform_curvature1*x(:,1)+transform_curvature2*x(:,2)) (x(:,2)+5).*sin(transform_curvature1*x(:,1)+transform_curvature2*x(:,2))];

data1 = f1(data_init);
data2 = f2(data_init);

data1_burst = zeros(size(data_burst_init));
data2_burst = zeros(size(data_burst_init));
for i=1:nsteps
    data1_burst(:, :, i) = f1(data_burst_init(:, :, i));
    data2_burst(:, :, i) = f2(data_burst_init(:, :, i));
end

%%

make_fig;
scatter(data_init(:,1),data_init(:,2),50,t,'.')
xlabel('x')
ylabel('y')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal
print('data_fastslow', fmt, res);

make_fig;
scatter(data1(:,1),data1(:,2),50,t,'.')
xlabel('x')
ylabel('z')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal
print('data_rescaled', fmt, res);

make_fig;
scatter(data2(:,1),data2(:,2),50,t,'.')
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal
print('data_nonlinear', fmt, res);

%% DMAPS

eps_scale = 100;

W = squareform(pdist(data_init)).^2;
eps = median(W(:));
[V_init_dmaps, D_init_dmaps] = dmaps(W, eps, 10);

make_fig;
scatter(data_init(:,1),data_init(:,2),50,V_init_dmaps(:,2),'.')
xlabel('x')
ylabel('y')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal
print('data_fastslow_dmaps', fmt, res);

W = squareform(pdist(data1)).^2;
eps_dmaps = median(W(:))/eps_scale;
[V1_dmaps, D1_dmaps] = dmaps(W, eps_dmaps, 10);

make_fig;
scatter(data1(:,1),data1(:,2),50,V1_dmaps(:,2),'.')
xlabel('x')
ylabel('z')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal
print('data_rescaled_dmaps', fmt, res);


W = squareform(pdist(data2)).^2;
eps = median(W(:));
[V2_dmaps, D2_dmaps] = dmaps(W, eps, 10);

make_fig;
scatter(data2(:,1),data2(:,2),50,V2_dmaps(:,2),'.')
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal
% print('data1_dmaps', fmt, res);

%%

idx = zeros(3, 1);
[~, idx(1)] = min(data_init(:,1));
[~, idx(2)] = max(data_init(:,1));
idx(3) = 1000;

make_fig;
plot(data_init(:,1),data_init(:,2),'.')
hold on
for i=1:3
    plot(data_burst_init(:, 1, idx(i)), data_burst_init(:,2,idx(i)),'.r')
end
xlabel('x')
ylabel('y')
axis equal 
print('data_fastslow_withbursts', fmt, res);

make_fig;
plot(data1(:,1),data1(:,2),'.')
hold on
for i=1:3
    plot(data1_burst(:, 1, idx(i)), data1_burst(:,2,idx(i)),'.r')
end
xlabel('x')
ylabel('z')
axis equal 
print('data_rescaled_withbursts', fmt, res);


make_fig;
plot(data2(:,1),data2(:,2),'.')
hold on
for i=1:3
    plot(data2_burst(:, 1, idx(i)), data2_burst(:,2,idx(i)),'.r')
end
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
axis equal 
print('data_nonlinear_withbursts', fmt, res);


%%
[inv_c_init, new_data_init, ranks_init] = covariances2(data_burst_init, dim);
[inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);
[inv_c2, new_data2, ranks2] = covariances2(data2_burst, dim);

%%

% eps_niv = eps_dmaps/(epsilon)^2;
eps_niv = eps_dmaps/(epsilon/200);
[V_init_niv, D_init_niv, eps, Dis_data_init] = NIV_return_dist(data_init, inv_c_init, eps_niv, 10, 0);


make_fig;
scatter(data_init(:,1),data_init(:,2),50,V_init_niv(:,2),'.')
xlabel('x')
ylabel('y')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal
print('data_fastslow_niv', fmt, res);

make_fig;
plot(V_init_niv(:,2), V1_dmaps(:,2),'.')
xlabel('NIV from fast-slow system')
ylabel('DMAPS from rescaled system')
axis equal
print('data_fastslow_corr', fmt, res);

%%


[V2_niv, D2_niv, eps, Dis_data2] = NIV_return_dist(data2, inv_c2, eps_niv, 10, 0);


make_fig;
scatter(data2(:,1),data2(:,2),50,V2_niv(:,2),'.')
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal
print('data_nonlinear_niv', fmt, res);

make_fig;
plot(V2_niv(:,2), V1_dmaps(:,2),'.')
xlabel('NIV from nonlinear fast-slow system')
ylabel('DMAPS from rescaled system')
axis equal
print('data_nonlinear_corr', fmt, res);

%%

W = squareform(pdist(data1)).^2;

figure;
plot(data_init(:,1),data_init(:,2),'.')
hold on
ind = find(Dis_data_init(:,1000) < eps_niv);
plot(data_init(ind,1),data_init(ind,2),'.r')
axis equal 


figure;
plot(data1(:,1),data1(:,2),'.')
hold on
plot(data1(ind,1),data1(ind,2),'.r')
axis equal

figure;
plot(data2(:,1),data2(:,2),'.')
hold on
ind = find(Dis_data2(:,1000) < eps_niv);
plot(data2(ind,1),data2(ind,2),'.r')

figure;
plot(data1(:,1),data1(:,2),'.')
hold on
ind = find(W(:,1000) < eps_dmaps);
plot(data1(ind,1),data1(ind,2),'.r')
axis equal

