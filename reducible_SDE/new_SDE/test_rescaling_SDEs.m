clear all
close all


%%

a = 1;
epsilon = 0.01;

DriftFn = @(t, x) [a; -x(2)/epsilon];
DiffnFn = @(t, x) [1 0; 0 1/sqrt(epsilon)];

tmax = 5;
dt = 1e-3;
dt_step = 1e-3;

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

mean(data1(:,2))
var(data1(:,2))

%%

a = 1;
epsilon = 0.01;

DriftFn = @(t, x) [a; -x(2)/epsilon];
DiffnFn = @(t, x) [1 0; 0 1];

tmax = 5;
dt = 1e-3;
dt_step = 1e-3;

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
scatter(data1(:,1),data1(:,2)/sqrt(epsilon),50,t)
axis equal

mean(data1(:,2)/sqrt(epsilon))
var(data1(:,2)/sqrt(epsilon))