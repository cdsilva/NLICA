clear all
close all

%%

N = 1500; % number of points considered

% archemedian spiral

a = 1;

theta_vec = linspace(0, 4*pi, 100);
s = 0.5*a*(theta_vec.*sqrt(1+theta_vec.^2)+log(theta_vec+sqrt(1+theta_vec.^2)));

theta = interp1(s, theta_vec, rand(N, 1)*max(s));

h = 40;

neigs = 10;

%%

rng(321);

z = h*rand(N,1); % random heights
x = a * cos(theta) .* theta;
y = a * sin(theta) .* theta;


%%

data = [x y z]; % data

W = squareform(pdist(data)).^2;
eps = 5;

[V, D] = dmaps(W, eps, neigs);

eps_med_scale = 3;
res = compute_residuals_DMAPS(V, eps_med_scale);

%%

make_fig(2,2);
plot3(x,y,z,'.', 'markersize', 10)
view(-20, 80)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
print('swissroll1.eps','-depsc')

make_fig(3,2);
scatter3(x,y,z,50,sign(corr(V(:,2),theta))*V(:,2),'.')
view(-20, 80)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_1');
print('swissroll1_color1.eps','-depsc')

make_fig(3,2);
scatter3(x,y,z,50,sign(corr(V(:,3),z))*V(:,3),'.')
view(-20, 80)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_2');
print('swissroll1_color2.eps','-depsc')

make_fig(3,2);
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'ylim', [0.9 1])
xlabel('k')
ylabel('\mu_k')
axis square
print('swissroll1_evals.eps','-depsc')

%%

z = z / 2;
data = [x y z]; % data

W = squareform(pdist(data)).^2;
eps = 4;

[V, D] = dmaps(W, eps, neigs);

eps_med_scale = 3;
res = compute_residuals_DMAPS(V, eps_med_scale);

%%

make_fig(2,2);
plot3(x,y,z,'.', 'markersize', 10)
view(-20, 80)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
% print('swissroll2.eps','-depsc')

make_fig(3,2);
scatter3(x,y,z,50,sign(corr(V(:,2),theta))*V(:,2),'.')
view(-20, 80)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_1');
% print('swissroll2_color1.eps','-depsc')

make_fig(3,2);
scatter3(x,y,z,50,sign(corr(V(:,6),z))*V(:,6),'.')
view(-20, 80)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_5');
% print('swissroll2_color2.eps','-depsc')

make_fig(3,2);
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'ylim', [0.9 1])
xlabel('k')
ylabel('\mu_k')
axis square
% print('swissroll2_evals.eps','-depsc')

