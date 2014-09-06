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

make_fig(3, 3);
scatter(data1(:,1),data1(:,2),50,t, '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
h = colorbar;
xlabel(h,'t');
saveas(gcf, 'data_init_nonlinear.eps', 'epsc');

%%

W = squareform(pdist(data1)).^2;
sigma = median(W(:))*100;
[V, D] = dmaps(W, sigma, 10);

make_fig(3, 3);
scatter(data1(:,1),data1(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
h = colorbar;
xlabel(h,'\phi_1');
saveas(gcf, 'data_nonlinear_DMAPS.eps', 'epsc');

%%

dim = 2;
[c, inv_c, ranks] = covariances(data1_burst, dim, dt_burst);

[V, D, eps, Dis] = NIV_return_dist(data1, inv_c, 0.01, 10, 0);

make_fig(3, 3);
scatter(data1(:,1),data1(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
h = colorbar;
xlabel(h,'\phi_1');
saveas(gcf, 'data_nonlinear_NIV.eps', 'epsc');

%%

dt_tmp = 10.^(-9:-2);
dt_tmp2 = logspace(-9, -2, 100);

make_fig(3, 3);
loglog(dt_tmp2, (1+4/epsilon+1/epsilon)*ones(size(dt_tmp2)), '-b')
hold on
loglog(dt_tmp2, 8*sqrt(dt_tmp2)/epsilon^(3/2), '-r')
loglog(dt_tmp2, (1+4/epsilon+1/epsilon)+8*sqrt(dt_tmp2)/epsilon^(3/2), '-k')
axis tight
axis square
xlabel('\tau - t')
% ylabel('$\| C \|$', 'interpreter','latex')
h = legend('$\|C\|$', '$\| $e$_C\|$', '$\| \hat{C} \| = \| C\| + \|$e$_C \|$', 'location','southwest');
 set(h,'Interpreter','latex');
 set(h,'fontsize', 6);
saveas(gcf, 'C_dt_analytical_nonlinear.eps', 'epsc');

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
make_fig(3, 3);
loglog(dt_tmp, mean(norm_c, 2),'.')
% axis([10^-6 10^-2 10^0 10^5])
% loglog(dt_tmp2, (1+4/epsilon+1/epsilon)+8*sqrt(dt_tmp2)/epsilon^(3/2))
% hold on
% loglog(dt_tmp2, (1+4/epsilon+1/epsilon)-8*sqrt(dt_tmp2)/epsilon^(3/2))
% for j=1:ntest_points
%     loglog(dt_tmp, norm_c(:,j),'.r')
%     hold on
% end
% axis tight
axis square
xlabel('\tau - t')
ylabel('$\| \hat{C} \|$', 'interpreter','latex')
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

make_fig(3, 3);
loglog(dy, (1+4+epsilon)*dy.^2,'-b')
hold on
loglog(dy, -(epsilon-1)*dy.^4,'-r')
loglog(dy, (1+4+epsilon)*dy.^2-(epsilon-1)*dy.^4,'-k')
axis tight
axis square
xlabel('$\|\vec{y}_2 - \vec{y}_1\|_2$', 'interpreter','latex')
% ylabel('$\|\vec{y}_2 - \vec{y}_1\|_M^2$', 'interpreter','latex')
h = legend('$\|\vec{y}_2 - \vec{y}_1\|_M^2$', '$\| $e$_M\|$', sprintf('%s\n%s','$\| \vec{z}_2 - \vec{z}_1 \|_2^2 = $','$ \|\vec{y}_2 - \vec{y}_1\|_M^2 + \|$e$_M \|$'), 'location','southeast');
 set(h,'Interpreter','latex');
 set(h,'fontsize', 4);
saveas(gcf, 'dist_dy_analytical_nonlinear.eps', 'epsc');


%%
W = squareform(pdist(data1));
minW = min(W(W > 0));
maxW = max(W(:));

dy2 = logspace(log10(minW), log10(maxW), 51);
dy2 = dy2(1:end-1);
Dis_avg = zeros(size(dy2));
for i=1:length(dy2)
    if i == length(dy2)
        Dis_avg(i) = mean(Dis(W >= dy2(i)));
    else      
        Dis_avg(i) = mean(Dis(W >= dy2(i) & W < dy2(i+1)));
    end
end

make_fig(3, 3);
rand_idx = randperm(nsteps^2, 2000);
loglog(dy2, Dis_avg, '.')
% hold on
% loglog(dy, (sqrt(epsilon)*(1+4+epsilon)*dy.^2+(epsilon-1)*dy.^4))
% hold on
% loglog(dy, (sqrt(epsilon)*(1+4+epsilon)*dy.^2-(epsilon-1)*dy.^4))
axis tight
axis square
xlabel('$\| \vec{y}_2 - \vec{y}_1 \|_2$', 'interpreter','latex')
ylabel('$\| \vec{y}_2 - \vec{y}_1 \|_M^2$', 'interpreter','latex')
saveas(gcf, 'dist_dy_nonlinear.eps', 'epsc');

%%

dt_burst = 1e-9;
eps_NIV = 10^-2;

[data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon);

data1 = f(data_init);
data1_burst = zeros(size(data_burst_init));
for i=1:nsteps
    data1_burst(:,:,i) = f(data_burst_init(:,:,i));
end

[c, inv_c, ranks] = covariances(data1_burst, dim, dt_burst);

[V, D, eps, Dis] = NIV_return_dist(data1, inv_c, eps_NIV, 10, 0);

make_fig(3, 3);
scatter(data1(:,1),data1(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
h = colorbar;
xlabel(h,'\phi_1');
saveas(gcf, 'data_nonlinear_NIV_good.eps', 'epsc');

%%
eps_NIV = 10^1;

[V, D, eps, Dis] = NIV_return_dist(data1, inv_c, eps_NIV, 10, 0);

make_fig(3, 3);
scatter(data1(:,1),data1(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
h = colorbar;
xlabel(h,'\phi_1');
saveas(gcf, 'data_nonlinear_NIV_bad_kernel.eps', 'epsc');


%%

dt_burst = 1e-4;
eps_NIV = 1e-2;

[data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon, data_init);

data1 = f(data_init);
data1_burst = zeros(size(data_burst_init));
for i=1:nsteps
    data1_burst(:,:,i) = f(data_burst_init(:,:,i));
end

[c, inv_c, ranks] = covariances(data1_burst, dim, dt_burst);

[V, D, eps, Dis] = NIV_return_dist(data1, inv_c, eps_NIV, 10, 0);


make_fig(3, 3);
scatter(data1(:,1),data1(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
h = colorbar;
xlabel(h,'\phi_1');
saveas(gcf, 'data_nonlinear_NIV_bad_dt.eps', 'epsc');


