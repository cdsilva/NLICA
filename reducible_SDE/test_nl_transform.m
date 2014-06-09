clear all
close all


%%

a = 1;
epsilon = 0.01;

DriftFn = @(t, x) [a; -x(2)/epsilon];
DiffnFn = @(t, x) [1 0; 0 1/sqrt(epsilon)];

tmax = 4;
dt = 0.5e-5;
dt_step = 0.5e-5;

nsteps = ceil(tmax / dt);
if mod(nsteps, 1) ~= 0
    disp('dt is not a multiple of tmax');
    return;
end

nsteps_per_step = dt / dt_step;
if mod(nsteps_per_step, 1) ~= 0
    disp('Step dt is not a multiple of dt');
    return;
end

data_start = [0; 0];

%%

rng(123);

SDE = sde(DriftFn, DiffnFn, 'StartState', data_start);
[data1, t, Z] = SDE.simulate(nsteps, 'DeltaTime', dt, 'NSTEPS', nsteps_per_step);

%%

figure;
scatter(data1(:,1),data1(:,2),50,t)
axis equal 

%%


% f = @(x) [(x(:,2)+5).*cos(x(:,1)+transform_curvature*x(:,2)) (x(:,2)+5).*sin(x(:,1)+transform_curvature*x(:,2))];

% R = max(abs(data1(:,2)))+0.01;
%f = @(x) [x(:,1)+sqrt(R^2-(x(:,2)).^2) x(:,2)];

freq = 0.5;
scale = 2;
f = @(x) [x(:,1)+scale*sin(x(:,2)/freq) x(:,2)];

data2 = f(data1);

figure;
scatter(data2(:,1),data2(:,2),50,t)
axis equal 

%% DMAPS

D = 2;
knn = 1000;
stride = ceil(size(data1, 1) / 3000);

%%
W = squareform(pdist(data1(1:stride:end,:))).^2;
eps = median(W(:));

[V1_dmaps, D1_dmaps] = dmaps(W, eps, 10);

figure;
scatter(data1(1:stride:end,1),data1(1:stride:end,2),50,V1_dmaps(:,2))
axis equal


%%
W = squareform(pdist(data2(1:stride:end,:))).^2;
eps = median(W(:));

[V2_dmaps, D2_dmaps] = dmaps(W, eps, 10);

figure;
scatter(data2(1:stride:end,1),data2(1:stride:end,2),50,V2_dmaps(:,2))
axis equal

figure;
scatter(data1(1:stride:end,1),data1(1:stride:end,2),50,V2_dmaps(:,2))
axis equal


figure;
plot(V1_dmaps(:,2), V2_dmaps(:,2),'.')


corr(V1_dmaps(:,2), V2_dmaps(:,2))


%%

[inv_c1, new_data1, ranks1] = covariances_diff(data1, knn, D, stride);
[inv_c2, new_data2, ranks2] = covariances_diff(data2, knn, D, stride);

%%

[V1_niv, D1_niv] = NIV(new_data1, inv_c1, 0, 10, 0);


figure;
scatter(new_data1(:,1),new_data1(:,2),50,V1_niv(:,2))
axis equal

%%

[V2_niv, D2_niv] = NIV(new_data2, inv_c2, 0, 10, 0);

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
