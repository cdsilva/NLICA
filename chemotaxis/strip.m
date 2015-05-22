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

%%

n = 2000;

rng(123);

for Lx = [2 4 8]
    
    data = rand(n, 2);
    data(:,1) = data(:,1) * Lx;
    data(:,2) = data(:,2) * Ly;
    
    W = squareform(pdist(data)).^2;
    eps = 0.15^2;
    
    [V, D] = dmaps(W, 2*eps, 10);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
    [~, idx] = sort(res, 'descend');
    idx = sort(idx(1:2));
    
    sqrt(log(D(idx(2), idx(2)))/log(D(idx(1), idx(1))))
    
    m = 0:10;
    m1 = reshape(repmat(m, length(m), 1), [], 1);
    m2 = reshape(repmat(m', 1, length(m)), [], 1);
    analytic_D = (m1/Lx).^2 + (m2/Ly).^2;
    [analytic_D, I] = sort(analytic_D);
    m1 = m1(I);
    m2 = m2(I);
    
    make_fig(3,2);
    scatter(data(:,1),data(:,2),25, 'b','.')
    axis square
    axis equal
    xlabel('z_1')
    ylabel('z_2')
    print(sprintf('strip_data_L%d.eps', Lx),'-depsc')
    
    
    make_fig(3,2);
    make_colored_bars(diag(D(2:end,2:end)), res(2:end))
    hold on
    plot(exp(-eps*pi^2*analytic_D(2:size(D,1))/2),'xr')
    xlabel('k')
    ylabel('\mu_k')
    axis square
    if Lx ~= 8
        colorbar off
    end
    print(sprintf('strip_spectrum_L%d.eps', Lx),'-depsc')
    
    
end

%%


figure;

n = 2000;
rng(321);

Lx = 4;
Ly = 1;

data = rand(n, 2);
data(:,1) = data(:,1) * Lx;
data(:,2) = data(:,2) * Ly;

plot_idx = 1;

subplot(2,4,plot_idx)
plot_idx = plot_idx + 1;
scatter(data(:,1),data(:,2),25, 'b','.')
axis square
axis equal
xlabel('z_1')
ylabel('z_2')

W = squareform(pdist(data)).^2;
eps = 0.15^2;
eps_med_scale = 3;

for alpha=[0 0.5 1]
    [V, D] = dmaps_weight(W, alpha, 2*eps, 10);
    res = compute_residuals_DMAPS(V, eps_med_scale);
    [~, idx] = sort(res, 'descend');
    idx = sort(idx(1:2));
    
    subplot(2,4,plot_idx)
    plot_idx = plot_idx + 1;
    make_colored_bars(diag(D(2:end,2:end)), res(2:end))
%     title(sprintf('\\alpha = %0.1f, L = %2.2f', alpha, sqrt(log(D(idx(2), idx(2)))/log(D(idx(1), idx(1))))))
    title(sprintf('\\alpha=%0.1f, L_x=%2.2f, L_y=%2.2f', alpha, pi/sqrt(log(D(idx(1), idx(1)))*(-2/eps)), pi/sqrt(log(D(idx(2), idx(2)))*(-2/eps))))
end



% data_tmp = (randn(10*n, 1)/6 + 0.5) * Lx;
data_tmp = (randn(10*n, 1)/4 + 0.5) * Lx;
idx = find(data_tmp > 0 & data_tmp < Lx);
data(:,1) = data_tmp(idx(1:n));

subplot(2,4,plot_idx)
plot_idx = plot_idx + 1;scatter(data(:,1),data(:,2),25, 'b','.')
axis square
axis equal
xlabel('z_1')
ylabel('z_2')

W = squareform(pdist(data)).^2;
eps = 0.15^2;
eps_med_scale = 3;

for alpha=[0 0.5 1]
    [V, D] = dmaps_weight(W, alpha, 2*eps, 10);
    res = compute_residuals_DMAPS(V, eps_med_scale);
    [~, idx] = sort(res, 'descend');
    idx = sort(idx(1:2));
    
    subplot(2,4,plot_idx)
    plot_idx = plot_idx + 1;
    make_colored_bars(diag(D(2:end,2:end)), res(2:end))
%     title(sprintf('\\alpha = %0.1f, L = %2.2f', alpha, sqrt(log(D(idx(2), idx(2)))/log(D(idx(1), idx(1))))))
    title(sprintf('\\alpha=%0.1f, L_x=%2.2f, L_y=%2.2f', alpha, pi/sqrt(log(D(idx(1), idx(1)))*(-2/eps)), pi/sqrt(log(D(idx(2), idx(2)))*(-2/eps))))

end

