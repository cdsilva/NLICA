clear all
close all

%%

N = 1500; % number of points considered
t = rand(1,N);
t = sort(4*pi*sqrt(t))'; 

h = 15*pi;

z = h*rand(N,1); % random heights
x = (t+.1).*cos(t);
y = (t+.1).*sin(t);
data = [x,y,z]; % data

W = squareform(pdist(data)).^2;
eps = 5;

neigs = 10;

[V, D] = dmaps(W, eps, neigs);

eps_med_scale = 3;
res = compute_residuals_DMAPS(V, eps_med_scale);

%%


figure;
plot3(x,y,z,'.', 'markersize', 10)
view(-20, 70)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on

figure;
scatter3(x,y,z,50,sign(corr(V(:,2),t))*V(:,2),'.')
view(-20, 75)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_1');

figure;
scatter3(x,y,z,50,sign(corr(V(:,3),z))*V(:,3),'.')
view(-20, 75)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_2');

figure;
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'ylim', [0.9 1])
xlabel('k')
ylabel('\mu_k')
axis square

%%

h = 6*pi;

z = h*rand(N,1); % random heights

data = [x,y,z]; % data

W = squareform(pdist(data)).^2;
eps = 2;

neigs = 10;

[V, D] = dmaps(W, eps, neigs);

eps_med_scale = 3;
res = compute_residuals_DMAPS(V, eps_med_scale);

%%


figure;
plot3(x,y,z,'.', 'markersize', 10)
view(-20, 70)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on

figure;
scatter3(x,y,z,50,sign(corr(V(:,2),t))*V(:,2),'.')
view(-20, 75)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_1');

figure;
scatter3(x,y,z,50,sign(corr(V(:,5),z))*V(:,5),'.')
view(-20, 75)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_2');

figure;
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'ylim', [0.9 1])
xlabel('k')
ylabel('\mu_k')
axis square



