clear all
close all

%%
rng(123);

n = 2000;

theta1 = 2*pi*rand(n, 1);
theta2 = 2*pi*rand(n, 1);

neigs = 20;

%%

r1 = 3;
r2 = 1;

x = (r1 + r2 * cos(theta2)) .* cos(theta1);
y = (r1 + r2 * cos(theta2)) .* sin(theta1);
z = r2 * sin(theta2);

data = [x y z];
W = squareform(pdist(data));
eps = median(W(:))/3;

[V, D] = dmaps(W.^2, eps.^2, neigs);

eps_med_scale = 3;
res = compute_residuals_DMAPS(V, eps_med_scale);


%
make_fig(3,2);
plot3(x,y,z,'.', 'markersize', 8)
% view(-20, 70)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
print('torus1.eps','-depsc')

make_fig(3,2);
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'xtick', (6:5:size(D,1))-1)
xlabel('k')
ylabel('\mu_k')
axis square
colorbar off
print('torus1_evals.eps','-depsc')

%%

r1 = 5;
r2 = 1;

x = (r1 + r2 * cos(theta2)) .* cos(theta1);
y = (r1 + r2 * cos(theta2)) .* sin(theta1);
z = r2 * sin(theta2);

data = [x y z];
W = squareform(pdist(data));
eps = median(W(:))/3;

[V, D] = dmaps(W.^2, eps.^2, neigs);

eps_med_scale = 3;
res = compute_residuals_DMAPS(V, eps_med_scale);

make_fig(3,2);
plot3(x,y,z,'.', 'markersize', 8)
% view(-20, 75)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
print('torus2.eps','-depsc')

make_fig(3,2);
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'xtick', (6:5:size(D,1))-1)
xlabel('k')
ylabel('\mu_k')
axis square
colorbar off
print('torus2_evals.eps','-depsc')

%%
r1 = 10;
r2 = 1;

x = (r1 + r2 * cos(theta2)) .* cos(theta1);
y = (r1 + r2 * cos(theta2)) .* sin(theta1);
z = r2 * sin(theta2);

data = [x y z];
W = squareform(pdist(data));
eps = median(W(:))/3;

[V, D] = dmaps(W.^2, eps.^2, neigs);

eps_med_scale = 3;
res = compute_residuals_DMAPS(V, eps_med_scale);

make_fig(3,2);
plot3(x,y,z,'.', 'markersize', 8)
% view(-20, 70)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
print('torus3.eps','-depsc')

make_fig(3,2);
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'xtick', (6:5:size(D,1))-1)
xlabel('k')
ylabel('\mu_k')
axis square
print('torus3_evals.eps','-depsc')


