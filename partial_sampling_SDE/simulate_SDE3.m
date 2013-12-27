function [res, x, y, energy_profile] = simulate_SDE3(nperiods, dt, D)

centers = [0.25 0.25;
   0.75 0.75;
   0.25 0.75];
% centers = [0.75 0.15;
%     0.25 0.5;
%     0.75 0.85];
scale = 0.05;

SDE = sde(@(t, X) DriftRate(t, X, centers, scale), @(t, X) DiffusionRate(t, X, D),'startstate',0.5*ones(2,1));

res = SDE.simulate(nperiods, 'DeltaTime', dt);
res = res(1:end-1,:);

x = 0:0.01:1;
y = 0:0.01:1;
x = x(1:end-1);
y = y(1:end-1);
[X, Y] = meshgrid(x, y);
energy_profile = energy(X, Y, centers, scale);

function Y = DriftRate(t, X, centers, scale)

[m, n] = size(centers);

x = X(1);
y = X(2);

Y = zeros(size(X));
norm_factor = 0;
for i=1:m
    Y(1) = Y(1) + exp(-((x-centers(i,1)).^2 + (y-centers(i,2)).^2)/scale)*(-2*(x-centers(i,1))/scale);
    Y(2) = Y(2) + exp(-((x-centers(i,1)).^2 + (y-centers(i,2)).^2)/scale)*(-2*(y-centers(i,2))/scale);
    norm_factor = norm_factor + exp(-((x-centers(i,1)).^2 + (y - centers(i,2)).^2)/scale);
end
Y = Y/norm_factor;

function Y = DiffusionRate(t, X, D)

Y = D*eye(length(X));

function f = energy(X,Y, centers, scale)

[m, n] = size(centers);

f = 0;
for i=1:m
    f = f + exp(-((X-centers(i,1)).^2 + (Y - centers(i,2)).^2)/scale);
end

f = -log(f);