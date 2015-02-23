clear all
close all

%%
rng(123);

n = 1000;

theta1 = 2*pi*rand(n, 1);
theta2 = 2*pi*rand(n, 1);

neigs = 20;

%%

% r1 = 3;
% r2 = 1;
% 
% x = (r1 + r2 * cos(theta2)) .* cos(theta1);
% y = (r1 + r2 * cos(theta2)) .* sin(theta1);
% z = r2 * sin(theta2);
% 
% data = [x y z];
% W = squareform(pdist(data));
% eps = median(W(:))/3;
% 
% [V, D] = dmaps(W.^2, eps.^2, neigs);
% 
% eps_med_scale = 3;
% res = compute_residuals_DMAPS(V, eps_med_scale);
% 
% save('torus1.mat')
load('torus1.mat')

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

% make_fig(3,3);
% scatter3(x,y,z,50,sign(corr(V(:,2),t))*V(:,2),'.')
% view(-20, 75)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% grid on
% cbar = colorbar('peer',gca); 
% set(get(cbar,'xlabel'),'String','\phi_1');
% print('swissroll1_color1.eps','-depsc')
% 
% make_fig(3,3);
% scatter3(x,y,z,50,sign(corr(V(:,3),z))*V(:,3),'.')
% view(-20, 75)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% grid on
% cbar = colorbar('peer',gca); 
% set(get(cbar,'xlabel'),'String','\phi_2');
% print('swissroll1_color2.eps','-depsc')

make_fig(3,2);
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'xtick', (6:5:size(D,1))-1)
xlabel('k')
ylabel('\mu_k')
axis square
print('torus1_evals.eps','-depsc')

%%

% r1 = 5;
% r2 = 1;
% 
% x = (r1 + r2 * cos(theta2)) .* cos(theta1);
% y = (r1 + r2 * cos(theta2)) .* sin(theta1);
% z = r2 * sin(theta2);
% 
% data = [x y z];
% W = squareform(pdist(data));
% eps = median(W(:))/3;
% 
% [V, D] = dmaps(W.^2, eps.^2, neigs);
% 
% eps_med_scale = 3;
% res = compute_residuals_DMAPS(V, eps_med_scale);
% 
% save('torus2.mat')
load('torus2.mat')

make_fig(3,2);
plot3(x,y,z,'.', 'markersize', 8)
% view(-20, 75)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
print('torus2.eps','-depsc')

% make_fig(3,3);
% scatter3(x,y,z,50,sign(corr(V(:,2),t))*V(:,2),'.')
% view(-20, 75)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% grid on
% cbar = colorbar('peer',gca); 
% set(get(cbar,'xlabel'),'String','\phi_1');
% print('swissroll2_color1.eps','-depsc')
% 
% make_fig(3,3);
% scatter3(x,y,z,50,sign(corr(V(:,6),z))*V(:,6),'.')
% view(-20, 75)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% grid on
% cbar = colorbar('peer',gca); 
% set(get(cbar,'xlabel'),'String','\phi_5');
% print('swissroll2_color2.eps','-depsc')


make_fig(3,2);
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'xtick', (6:5:size(D,1))-1)
xlabel('k')
ylabel('\mu_k')
axis square
print('torus2_evals.eps','-depsc')

%%
% r1 = 10;
% r2 = 1;
% 
% x = (r1 + r2 * cos(theta2)) .* cos(theta1);
% y = (r1 + r2 * cos(theta2)) .* sin(theta1);
% z = r2 * sin(theta2);
% 
% data = [x y z];
% W = squareform(pdist(data));
% eps = median(W(:))/3;
% 
% [V, D] = dmaps(W.^2, eps.^2, neigs);
% 
% eps_med_scale = 3;
% res = compute_residuals_DMAPS(V, eps_med_scale);
% 
% save('torus3.mat')
load('torus3.mat')

%
make_fig(3,2);
plot3(x,y,z,'.', 'markersize', 8)
% view(-20, 70)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
print('torus3.eps','-depsc')

% make_fig(3,3);
% scatter3(x,y,z,50,sign(corr(V(:,2),t))*V(:,2),'.')
% view(-20, 75)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% grid on
% cbar = colorbar('peer',gca); 
% set(get(cbar,'xlabel'),'String','\phi_1');
% print('swissroll1_color1.eps','-depsc')
% 
% make_fig(3,3);
% scatter3(x,y,z,50,sign(corr(V(:,3),z))*V(:,3),'.')
% view(-20, 75)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% grid on
% cbar = colorbar('peer',gca); 
% set(get(cbar,'xlabel'),'String','\phi_2');
% print('swissroll1_color2.eps','-depsc')

make_fig(3,2);
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'xtick', (6:5:size(D,1))-1)
xlabel('k')
ylabel('\mu_k')
axis square
print('torus3_evals.eps','-depsc')


