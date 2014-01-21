function [t, n_all] = SDE_simulation(n0, tmax, dt, nsteps_per_step)

parameters;

SDE = sde(@(t,x) DriftFn(t, x, k), @(t,x) DiffusionFn(t, x, k), 'StartState', n0);

nPeriods = ceil(tmax / dt);
[n_all, t, Z] = SDE.simulate(nPeriods, 'DeltaTime', dt, 'NSTEPS', nsteps_per_step);

% function Y = DriftFn(t, X, k)
% 
% r = rates(X, k);
% 
% Y = zeros(5, 1);
% 
% Y(1) = -r(1) + r(2);
% Y(2) = -r(1) + r(2);
% Y(3) = r(1) - r(2) - r(3);
% Y(4) = -r(3);
% Y(5) = r(3);
% 
% function Y = DiffusionFn(t, X, k)
% 
% r = rates(X, k);
% 
% Y = zeros(5, 3);
% 
% Y(1,:) = [-sqrt(r(1)) sqrt(r(2)) 0];
% Y(2,:) = [-sqrt(r(1)) sqrt(r(2)) 0];
% Y(3,:) = [sqrt(r(1)) -sqrt(r(2)) -sqrt(r(3))];
% Y(4,:) = [0 0 -sqrt(r(3))];
% Y(5,:) = [0 0 sqrt(r(3))];
% 
% function r = rates(X, k)
% 
% r = zeros(3, 1);
% r(1) = k(1)*X(1)* X(2);
% r(2) = k(2)*X(3);
% r(3) = k(3)*X(3)*X(4);

    
