clear all
close all

%%

epsilon = 1e-3;

dt_burst = 1e-9;
dt = 1e-4;

nsteps = 3000;
nsteps_burst = 50;

rng(123);
[data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon);

make_fig(6, 6);
scatter(data_init(:,1),data_init(:,2),50,t, '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
saveas(gcf, 'data_init.eps');


%%

W = squareform(pdist(data_init)).^2;
sigma = median(W(:));
[V, D] = dmaps(W, sigma, 10);

make_fig(6, 6);
scatter(data_init(:,1),data_init(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
saveas(gcf, 'data_linear_DMAPS.eps');

%%

dim = 2;
[c, inv_c, ranks] = covariances(data_burst_init, dim, dt_burst);

[V, D, eps, Dis] = NIV_return_dist(data_init, inv_c, 0, 10, 0);

make_fig(6, 6);
scatter(data_init(:,1),data_init(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
saveas(gcf, 'data_linear_NIV.eps');

%%

dt_tmp = 10.^(-8:-2);
dt_tmp2 = logspace(-8, -2, 100);

figure; 
loglog(dt_tmp2, (1+1/epsilon)*ones(size(dt_tmp2)), '-b')
hold on
loglog(dt_tmp2, dt_tmp2/epsilon^2, '-r')

figure; 
loglog(dt_tmp2, 1+1/epsilon+dt_tmp2/epsilon^2)
hold on
loglog(dt_tmp2, 1+1/epsilon-dt_tmp2/epsilon^2)

%%

ntest_points = 10;
data_init_tmp = data_init(1:10:10*ntest_points,:);
norm_c = zeros(length(dt_tmp), ntest_points);

for i=1:length(dt_tmp)
    [~, data_burst_init_tmp, ~] = simulate_SDE(nsteps, dt, nsteps_burst, dt_tmp(i), epsilon, data_init_tmp);
    [c, inv_c, ranks] = covariances(data_burst_init_tmp, dim, dt_tmp(i));
    for j=1:ntest_points
        norm_c(i, j) = norm(c(:,:,j));
    end
end

%%
make_fig(6, 6);
loglog(dt_tmp2, 1+1/epsilon+dt_tmp2/epsilon^2)
hold on
loglog(dt_tmp2, 1+1/epsilon-dt_tmp2/epsilon^2)
for j=1:ntest_points
    loglog(dt_tmp, norm_c(:,j),'.r')
    hold on
end
axis tight
axis square
xlabel('\tau - t')
ylabel('$\| C \|$', 'interpreter','latex')
saveas(gcf, 'C_dt_linear.eps', 'epsc');

%%

dt_burst = 1e-6;
[data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon, data_init);

[c, inv_c, ranks] = covariances(data_burst_init, dim, dt_burst);

[V, D, eps, Dis] = NIV_return_dist(data_init, inv_c, 0, 10, 0);

make_fig(6, 6);
scatter(data_init(:,1),data_init(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
saveas(gcf, 'data_linear_NIV_smalldt.eps', 'epsc');

%%

dt_burst = 1e-3;
[data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon, data_init);

[c, inv_c, ranks] = covariances(data_burst_init, dim, dt_burst);

[V, D, eps, Dis] = NIV_return_dist(data_init, inv_c, 0, 10, 0);

make_fig(6, 6);
scatter(data_init(:,1),data_init(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
saveas(gcf, 'data_linear_NIV_largedt.eps', 'epsc');





