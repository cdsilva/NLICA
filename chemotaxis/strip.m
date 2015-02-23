clear all
close all

%%

Lx = 4;
Ly = 1;


delta = 0.1;

[X, Y] = meshgrid(0:delta:Lx, 0:delta:Ly);

make_fig(3,2);
h = pcolor(X, Y, cos(pi*X/Lx));
set(h, 'edgecolor','none');
shading interp
axis equal
xlabel('z_1')
ylabel('z_2')
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','$\tilde{\phi}_{1,0}$', 'interpreter','latex');
print('strip_cnts1.eps','-depsc')

make_fig(3,2);
h = pcolor(X, Y, cos(pi*Y/Ly));
set(h, 'edgecolor','none');
shading interp
axis equal
xlabel('z_1')
ylabel('z_2')
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','$\tilde{\phi}_{0,1}$', 'interpreter','latex');
print('strip_cnts2.eps','-depsc')


%%

rng(321);

n = 1500;
data = rand(n, 2);
data(:,1) = data(:,1) * Lx;
data(:,2) = data(:,2) * Ly;

W = squareform(pdist(data)).^2;
% eps = 0.25^2;
eps = median(W(:));

[V, D] = dmaps(W, eps, 10);

if corr(V(:,2), data(:,1)) > 0
    V(:,2) = -V(:,2);
end
if corr(V(:,5), data(:,2)) > 0
    V(:,5) = -V(:,5);
end

for i=2:5
make_fig(3,2);
scatter(data(:,1),data(:,2),50, V(:,i)/max(V(:,i)),'.')
axis equal
xlabel('z_1')
ylabel('z_2')
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String',sprintf('\\phi_{%d}', i-1));
print(sprintf('strip_discrete%d.eps', i-1),'-depsc')
end

%%

rng(321);

data_tmp = (randn(10*n, 1)/6 + 0.5) * Lx;
idx = find(data_tmp > 0 & data_tmp < Lx);
data(:,1) = data_tmp(idx(1:n));

W = squareform(pdist(data)).^2;
% eps = 0.25^2;
eps = median(W(:));

[V, D] = dmaps(W, eps, 10);

if corr(V(:,2), data(:,1)) > 0
    V(:,2) = -V(:,2);
end
if corr(V(:,4), data(:,2)) > 0
    V(:,4) = -V(:,4);
end

make_fig(3,2);
scatter(data(:,1),data(:,2),50, V(:,2)/max(V(:,2)),'.')
axis equal
xlabel('z_1')
ylabel('z_2')
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_{1}');
print('strip_nonuniform1.eps','-depsc')

make_fig(3,2);
scatter(data(:,1),data(:,2),50, V(:,4)/max(V(:,4)),'.')
axis equal
xlabel('z_1')
ylabel('z_2')
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_{3}');
print('strip_nonuniform2.eps','-depsc')

