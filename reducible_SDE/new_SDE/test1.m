clear all
close all


%%

fmt = '-djpeg';
res = '-r600';

dim = 2;

a = 0.1;
epsilon = 0.01;

DriftFn = @(t, x) [a; -x(2)];
DiffnFn = @(t, x) [1 0; 0 1];

dt = 0.5;
dt_step = dt;
nsteps = 2000;

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

make_fig;
scatter(data_init(:,1),data_init(:,2),50,t,'.')
xlabel('x')
ylabel('y')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal 
% print('raw_data', fmt, res);

%%

nsteps_burst = 200;
dt_burst = dt/100;
dt_burst_step = dt_burst;

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

%%

W = squareform(pdist(data_init)).^2;
eps = median(W(:))/10;

[V_init_dmaps, D_init_dmaps] = dmaps(W, eps, 10);

make_fig;
scatter(data_init(:,1),data_init(:,2),50,V_init_dmaps(:,2),'.')
xlabel('x')
ylabel('y')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal 
% print('raw_data_dmaps', fmt, res);

%%

make_fig;
plot(data_init(:,1),data_init(:,2),'.')
hold on

idx = zeros(3, 1);
[~, idx(1)] = min(data_init(:,1));
[~, idx(2)] = max(data_init(:,1));
idx(3) = 1000;
for i=1:3
    plot(data_burst_init(:, 1, idx(i)), data_burst_init(:,2,idx(i)),'.r')
end
xlabel('x')
ylabel('y')
axis equal 
% print('raw_data_withbursts', fmt, res);

%% first function (scaling)

scale = 100;
f = @(x) [x(:,1) scale*x(:,2)];

data1 = f(data_init);

data1_burst = zeros(size(data_burst_init));
for i=1:nsteps
    data1_burst(:, :, i) = f(data_burst_init(:, :, i));
end

make_fig;
scatter(data1(:,1),data1(:,2),50,t,'.')
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal 
% print('data1', fmt, res);

%%

make_fig;
plot(data1(:,1),data1(:,2),'.')
hold on
for i=1:3
    plot(data1_burst(:, 1, idx(i)), data1_burst(:,2,idx(i)),'.r')
end
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
axis equal 
% print('data1_withbursts', fmt, res);

%% DMAPS

W = squareform(pdist(data1)).^2;
eps = median(W(:));

[V1_dmaps, D1_dmaps] = dmaps(W, eps, 10);

make_fig;
scatter(data1(:,1),data1(:,2),50,V1_dmaps(:,2),'.')
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal 
% print('data1_dmaps', fmt, res);

%%

[inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);

%%

[V1_niv, D1_niv, eps] = NIV(data1, inv_c1, 1e4, 10, 0);


make_fig;
scatter(data1(:,1),data1(:,2),50,V1_niv(:,2),'.')
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal 
% print('data1_niv', fmt, res);

%%

make_fig;
plot(V1_niv(:,2),V_init_dmaps(:,2),'.')
xlabel('NIV for data with large noise')
ylabel('DMAPS for data with small noise')
axis equal 
% print('data1_corr', fmt, res);

%% second function (nonlinear transform)

% transform_curvature1 = 0.05;
% transform_curvature2 = 0.25;
transform_curvature1 = 0.025;
transform_curvature2 = 0.25;
f2 = @(x) [(x(:,2)+5).*cos(transform_curvature1*x(:,1)+transform_curvature2*x(:,2)) (x(:,2)+5).*sin(transform_curvature1*x(:,1)+transform_curvature2*x(:,2))];

data2 = f2(data_init);
data2_burst = zeros(size(data_burst_init));
for i=1:nsteps
    data2_burst(:, :, i) = f2(data_burst_init(:, :, i));
end

make_fig;
scatter(data2(:,1),data2(:,2),50,t,'.')
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal 
% print('data2', fmt, res);

%%

make_fig;
plot(data2(:,1),data2(:,2),'.')
hold on
for i=1:3
    plot(data2_burst(:, 1, idx(i)), data2_burst(:,2,idx(i)),'.r')
end
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
axis equal 
% print('data2_withbursts', fmt, res);

%%

W = squareform(pdist(data2)).^2;
eps = median(W(:));

[V2_dmaps, D2_dmaps] = dmaps(W, eps, 10);

make_fig;
scatter(data2(:,1),data2(:,2),50,V2_dmaps(:,2))


make_fig;
plot(V2_dmaps(:,2),V_init_dmaps(:,2),'.')

%%

[inv_c2, new_data2, ranks2] = covariances2(data2_burst, dim);


%%

[V2_niv, D2_niv, eps, Dis_data2] = NIV_return_dist(data2, inv_c2, 1e4, 10, 0);

figure;
plot(data2(:,1),data2(:,2),'.')
hold on
ind = find(Dis_data2(:,1000) < eps);
plot(data2(ind,1),data2(ind,2),'.r')

make_fig;
scatter(data2(:,1),data2(:,2),50,V2_niv(:,2),'.')
xlabel('f_1(x,y)')
ylabel('f_2(x,y)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal 
% print('data2_niv', fmt, res);

%%

make_fig;
plot(V2_niv(:,2),V_init_dmaps(:,2),'.')
xlabel('NIV for data with curved noise')
ylabel('DMAPS for data with small noise')
axis equal 
% print('data2_corr', fmt, res);
