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

f = @(x) [x(:,1)+x(:,2).^2 x(:,2)];


data1 = f(data_init);
data1_burst = zeros(size(data_burst_init));
for i=1:nsteps
    data1_burst(:,:,i) = f(data_burst_init(:,:,i));
end

make_fig(6, 6);
scatter(data_init(:,1),data_init(:,2),50,t, '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')

make_fig(6, 6);
scatter(data1(:,1),data1(:,2),50,t, '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
saveas(gcf, 'data_init_nonlinear.eps');

%%

W = squareform(pdist(data1)).^2;
sigma = median(W(:))*100;
[V, D] = dmaps(W, sigma, 10);

make_fig(6, 6);
scatter(data1(:,1),data1(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
saveas(gcf, 'data_nonlinear_DMAPS.eps');

%%

dim = 2;
[c, inv_c, ranks] = covariances(data1_burst, dim, dt_burst);

[V, D, eps, Dis] = NIV_return_dist(data1, inv_c, 0.01, 10, 0);

make_fig(6, 6);
scatter(data1(:,1),data1(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
saveas(gcf, 'data_nonlinear_NIV.eps');

%%

dt_tmp = 10.^(-8:-2);
dt_tmp2 = logspace(-8, -2, 100);

figure; 
loglog(dt_tmp2, (1+1/epsilon)*ones(size(dt_tmp2)), '-b')
hold on
loglog(dt_tmp2, dt_tmp2/epsilon^2, '-r')

figure; 
loglog(dt_tmp2, (1+4/epsilon+1/epsilon)+8*sqrt(dt_tmp2)/epsilon^(3/2))
hold on
loglog(dt_tmp2, (1+4/epsilon+1/epsilon)-8*sqrt(dt_tmp2)/epsilon^(3/2))

%%

ntest_points = 10;
data_init_tmp = data_init(1:10:10*ntest_points,:);
data1_tmp = f(data_init_tmp);

norm_c = zeros(length(dt_tmp), ntest_points);

for i=1:length(dt_tmp)
    [~, data_burst_init_tmp, ~] = simulate_SDE(nsteps, dt, nsteps_burst, dt_tmp(i), epsilon, data_init_tmp);
    data1_burst_tmp = zeros(size(data_burst_init_tmp));
    for j=1:ntest_points
        data1_burst_tmp(:,:,j) = f(data_burst_init_tmp(:,:,j));
    end
    [c, inv_c, ranks] = covariances(data1_burst_tmp, dim, dt_tmp(i));
    for j=1:ntest_points
        norm_c(i, j) = norm(c(:,:,j));
    end
end

%%
make_fig(6, 6);
loglog(dt_tmp2, (1+4/epsilon+1/epsilon)+8*sqrt(dt_tmp2)/epsilon^(3/2))
hold on
loglog(dt_tmp2, (1+4/epsilon+1/epsilon)-8*sqrt(dt_tmp2)/epsilon^(3/2))
for j=1:ntest_points
    loglog(dt_tmp, norm_c(:,j),'.r')
    hold on
end
axis tight
axis square
xlabel('\tau - t')
ylabel('$\| C \|$', 'interpreter','latex')
saveas(gcf, 'C_dt_nonlinear.eps', 'epsc');


%%

dt_burst = 1e-9;
nsteps_burst = 25;
dt = 1e-4;

[data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon);

data1 = f(data_init);
data1_burst = zeros(size(data_burst_init));
for i=1:nsteps
    data1_burst(:,:,i) = f(data_burst_init(:,:,i));
end

[c, inv_c, ranks] = covariances(data1_burst, dim, dt_burst);

[V, D, eps, Dis] = NIV_return_dist(data1, inv_c, 0, 10, 0);


%%

dy = logspace(-4, 2, 100);

figure;
loglog(dy, epsilon*dy.^2+(epsilon-1)*dy.^4)
hold on
loglog(dy, epsilon*dy.^2-(epsilon-1)*dy.^4)

W = squareform(pdist(data1));

make_fig(6, 6);
rand_idx = randperm(nsteps^2, 2000);
loglog(W(rand_idx), Dis(rand_idx),'.r')
hold on
loglog(dy, 0.1*((1+4+epsilon)*dy.^2+(epsilon-1)*dy.^4))
hold on
loglog(dy, 0.1*((1+4+epsilon)*dy.^2-(epsilon-1)*dy.^4))
axis tight
axis square
xlabel('$\| Y_2 - Y_1 \|_2$', 'interpreter','latex')
ylabel('$\| Y_2 - Y_1 \|_M^2$', 'interpreter','latex')
saveas(gcf, 'dist_dy_nonlinear.eps', 'epsc');

return


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





