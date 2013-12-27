% function does two-welled potential in 1-D

function [res, x, energy_profile] = simulate_SDE2(nperiods, dt, D, initial_point)

centers = [0.25; 0.75];
scale = 0.05;

SDE = sde(@(t, X) DriftRate(t, X, centers, scale), @(t, X) DiffusionRate(t, X, D),'startstate',initial_point);

res = SDE.simulate(nperiods, 'DeltaTime', dt);
res = res(1:end-1,:);

x = 0:0.01:1;
x = x(1:end-1);
energy_profile = energy(x, centers, scale);

function Y = DriftRate(t, x, centers, scale)

m = length(centers);

Y = 0;
norm_factor = 0;
for i=1:m
    Y = Y - (2 * (x-centers(i)) / scale) * exp(-(x-centers(i)).^2/scale);
    norm_factor = norm_factor + exp(-(x-centers(i)).^2/scale);
end

Y = Y/norm_factor;

function Y = DiffusionRate(t, X, D)

Y = D;

function f = energy(X, centers, scale)

f = -log(p(X, centers, scale));

function y = p(x, centers, scale)

m = length(centers);

y = zeros(size(x));
for i=1:m
    y = y + exp(-(x-centers(i)).^2/scale);
end


    