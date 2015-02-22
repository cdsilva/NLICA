clear all
close all

%%

Lx = 3;
Ly = 1;


delta = 0.1;

[X, Y] = meshgrid(0:delta:Lx, 0:delta:Ly);

figure;
h = pcolor(X, Y, cos(pi*X/Lx));
set(h, 'edgecolor','none');
shading interp
axis equal

figure;
h = pcolor(X, Y, cos(pi*Y/Ly));
set(h, 'edgecolor','none');
shading interp
axis equal

%%

rng(321);

n = 1000;
data = rand(n, 2);
data(:,1) = data(:,1) * Lx;
data(:,2) = data(:,2) * Ly;

W = squareform(pdist(data)).^2;
eps = 0.5^2;

[V, D] = dmaps(W, eps, 10);

figure;
scatter(data(:,1),data(:,2),50, V(:,2),'.')
axis equal

figure;
scatter(data(:,1),data(:,2),50, V(:,4),'.')
axis equal