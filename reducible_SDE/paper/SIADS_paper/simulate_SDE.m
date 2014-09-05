function [data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon, data_init)

dim = 2;

a = 3;
DriftFn = @(t, x) [a; -x(2)/epsilon];
DiffnFn = @(t, x) [1 0; 0 1/sqrt(epsilon)];

% min_dt = 1e-5;
min_dt = epsilon^2;
step_tol = 1e-4;

if nargin < 6
    dt_step = min(dt, min_dt);
    nsteps_per_step = dt / dt_step;
    if abs(mod(nsteps_per_step, 1)) > step_tol
        disp('Step dt is not a multiple of dt');
        return;
    end
    nsteps_per_step = round(nsteps_per_step);

    data_start = [0; 0];

    SDE = sde(DriftFn, DiffnFn, 'StartState', data_start);
    [data_init, t, ~] = SDE.simulate(nsteps, 'DeltaTime', dt, 'NSTEPS', nsteps_per_step);
    data_init = data_init(2:end, :);
    t = t(2:end);
else
    nsteps = size(data_init, 1);
    t = (1:nsteps)';
end

dt_burst_step = min(dt_burst, min_dt);
nsteps_burst_per_step = round(dt_burst / dt_burst_step);
if abs(mod(nsteps_burst_per_step, 1)) > step_tol
    disp('Step dt_burst is not a multiple of dt_burst');
    return;
end
nsteps_burst_per_step = round(nsteps_burst_per_step);

data_burst_init = zeros(nsteps_burst, dim, nsteps);
for i=1:nsteps
    i
    SDE = sde(DriftFn, DiffnFn, 'StartState', data_init(i, :)');
    [data_tmp, ~, ~] = SDE.simulate(1, 'DeltaTime', dt_burst, 'NSTEPS', nsteps_burst_per_step, 'ntrials', nsteps_burst);
    data_burst_init(:, :, i) = squeeze(data_tmp(end, :, :))';
end

