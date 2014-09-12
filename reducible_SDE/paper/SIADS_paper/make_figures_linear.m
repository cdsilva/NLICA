clear all
close all

%%

epsilon = 1e-3;

dt_burst = 1e-6;
dt = 1e-4;

nsteps = 3000;
nsteps_burst = 50;

rng(123);
[data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon);

make_fig(3, 3);
scatter(data_init(:,1),data_init(:,2),50,t, '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
h = colorbar;
xlabel(h,'t');
saveas(gcf, 'data_init.eps', 'epsc');

%%

W = squareform(pdist(data_init)).^2;
sigma = median(W(:));
[V, D] = dmaps(W, sigma, 10);

V(:,2) = sign(corr(V(:,2), data_init(:,2)))*V(:,2);

make_fig(3, 3);
scatter(data_init(:,1),data_init(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
h = colorbar;
xlabel(h,'\phi_1');
saveas(gcf, 'data_linear_DMAPS.eps','epsc');

%%

dim = 2;
[c, inv_c, ranks] = covariances(data_burst_init, dim, dt_burst);

[V, D, eps, Dis] = NIV_return_dist(data_init, inv_c, 0, 10, 0);
V(:,2) = sign(corr(V(:,2), data_init(:,1)))*V(:,2);

make_fig(3, 3);
scatter(data_init(:,1),data_init(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('x_1')
ylabel('x_2')
h = colorbar;
xlabel(h,'\phi_1');
saveas(gcf, 'data_linear_NIV.eps','epsc');

%%

dt_tmp = 10.^(-6:-2);
dt_tmp2 = logspace(-7, -2, 100);
dt_burst_samplepoints = [1e-6 1e-5 1e-3];

make_fig(3, 3);
loglog(dt_tmp2, (1+1/epsilon)*ones(size(dt_tmp2)), '-b')
hold on
% loglog(dt_tmp2, 2*dt_tmp2/epsilon^2, '-r')
loglog(dt_tmp2, dt_tmp2/(epsilon^2), '-r')
loglog(dt_tmp2, 1+1/epsilon+dt_tmp2/(epsilon^2), '-k')
for i=1:length(dt_burst_samplepoints)
    loglog(dt_burst_samplepoints(i)*ones(size(dt_tmp2)), dt_tmp2/epsilon^2, 'linestyle', ':', 'color', 0.5*ones(1,3))
end
axis([10^-7 10^-2 10^0 10^4])
axis square
xlabel('\delta t')
ylabel('Analytical covariance contributions')
h = legend('$\|C\|$', '$\| $e$_C\|$', '$\| C\| + \|$e$_C \|$', 'location','southwest');
set(h,'Interpreter','latex');
set(h,'fontsize', 6);
saveas(gcf, 'C_dt_analytical_linear.eps', 'epsc');

return

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
make_fig(3.5, 3);
loglog(dt_tmp, mean(norm_c, 2),'.')
hold on
for i=1:length(dt_burst_samplepoints)
    loglog(dt_burst_samplepoints(i)*ones(size(dt_tmp2)), dt_tmp2/epsilon^2, 'linestyle', ':', 'color', 0.5*ones(1,3))
end
axis([10^-7 10^-2 10^0 10^4])
% loglog(dt_tmp2, 1+1/epsilon+dt_tmp2/epsilon^2)
% hold on
% loglog(dt_tmp2, 1+1/epsilon-dt_tmp2/epsilon^2)
% for j=1:ntest_points
%     loglog(dt_tmp, norm_c(:,j),'.r')
%     hold on
% end
% % axis tight
% axis([10^-8 10^-1 -5000 5000])
axis square
xlabel('\delta t')
ylabel('$\| \hat{C} \|$', 'interpreter','latex')
saveas(gcf, 'C_dt_linear.eps', 'epsc');

%%

for i=1:length(dt_burst_samplepoints)
    dt_burst = dt_burst_samplepoints(i);
    [data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon, data_init);
    
    [c, inv_c, ranks] = covariances(data_burst_init, dim, dt_burst);
    
    [V, D, eps, Dis] = NIV_return_dist(data_init, inv_c, 0, 25, 0);
    
    [~, idx] = max(abs(corr(V, data_init(:,2))));
    
    V(:,2) = sign(corr(V(:,2),data_init(:,1)))*V(:,2);
    V(:,idx) = sign(corr(V(:,idx),data_init(:,2)))*V(:,idx);
    
    make_fig(2, 2);
    plot(data_init(:,1), data_init(:,2), '.', 'color', 0.5*ones(1,3))
    hold on
    plot(data_burst_init(:,1,1000), data_burst_init(:,2,1000),'.r')
    axis tight
    axis square
    xlabel('x_1')
    ylabel('x_2')
    saveas(gcf, sprintf('data_linear_burst%d.eps', i), 'epsc');
    
    make_fig(2, 2);
    plot(data_init(:,1), V(:, 2), '.')
    axis tight
    axis square
    xlabel('x_1')
    ylabel('\phi_1')
    saveas(gcf, sprintf('data_linear_slow%d.eps', i), 'epsc');
    
    make_fig(2, 2);
    plot(data_init(:,2), V(:, idx), '.')
    axis tight
    axis square
    xlabel('x_2')
    ylabel(sprintf('\\phi_{%d}', idx-1))
    saveas(gcf, sprintf('data_linear_fast%d.eps', i), 'epsc');
    
    make_fig(2, 2);
    plot(abs(diag(D(2:end, 2:end))),'.')
    hold on
    plot(1, abs(D(2,2)),'or')
    plot(idx-1, abs(D(idx, idx)), 'or')
    axis tight
    axis square
    axis([0 25 0 1])
    xlabel('k')
    ylabel('\lambda_k')
    saveas(gcf, sprintf('data_linear_evals%d.eps', i), 'epsc');
    
end







