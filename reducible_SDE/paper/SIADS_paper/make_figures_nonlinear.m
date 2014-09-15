clear all
close all

%%

epsilon = 1e-3;

dt_burst = 1e-9;
dt = 1e-4;

nsteps = 3000;
nsteps_burst = 50;

f = @(x) [x(:,1)+x(:,2).^2 x(:,2)];
dim = 2;

kernel_samplepoints = [1e-2 1e1];
dt_burst_samplepoints = [1e-7 1e-3];

%%
rng(123);
[data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon);

data1 = f(data_init);
data1_burst = zeros(size(data_burst_init));
for i=1:nsteps
    data1_burst(:,:,i) = f(data_burst_init(:,:,i));
end

%%
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
xlabel('y_1')
ylabel('y_2')
h = colorbar;
xlabel(h,'t');
saveas(gcf, 'data_init_nonlinear.eps', 'epsc');

%%

W = squareform(pdist(data1)).^2;
% sigma = median(W(:))*100;
sigma = 1e1;
[V, D] = dmaps(W, sigma, 10);
V(:,2) = sign(corr(V(:,2), data_init(:,2)))*V(:,2);

make_fig(3, 3);
scatter(data1(:,1),data1(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('y_1')
ylabel('y_2')
h = colorbar;
xlabel(h,'\phi_1');
saveas(gcf, 'data_nonlinear_DMAPS.eps', 'epsc');

%%

[c, inv_c, ranks] = covariances(data1_burst, dim, dt_burst);

[V, D, eps, Dis] = NIV_return_dist(data1, inv_c, 1e-2, 10, 0);
V(:,2) = sign(corr(V(:,2), data_init(:,1)))*V(:,2);

make_fig(3, 3);
scatter(data1(:,1),data1(:,2),50,V(:,2), '.')
axis tight
axis square
xlabel('y_1')
ylabel('y_2')
h = colorbar;
xlabel(h,'\phi_1');
saveas(gcf, 'data_nonlinear_NIV.eps', 'epsc');

%%
% Dis_analytical = zeros(size(Dis));
% for i=1:200
%     for j=1:i-1
%         invc_tmp = [1 -2*data_init(i,2); -2*data_init(i,2) epsilon+4*data_init(i,2)^2];
%         Dis_analytical(i,j) = (data1(i,:)-data1(j,:))*invc_tmp*(data1(i,:)-data1(j,:))';
%     end
% end
% Dis_analytical = 0.5*(Dis_analytical + Dis_analytical');
% figure; loglog(Dis(1:200,1), Dis_analytical(1:200,1),'.')
% hold on
% loglog(Dis(1:200,1), Dis(1:200,1),'.r')

%%

dy = logspace(-4, 1.5, 100);

make_fig(3, 3.5);
loglog(dy, 1.5*dy.^2,'-b')
hold on
loglog(dy, 1*dy.^4,'-r')
loglog(dy, 1.5*dy.^2+1*dy.^4,'-k')
for i=1:length(kernel_samplepoints)
    loglog(dy, kernel_samplepoints(i)*ones(size(dy)), 'linestyle', ':', 'color', 0.5*ones(1,3))
end
axis([1e-2 1e1 1e-6 1e2])
axis square
xlabel('$\|\vec{y}_2 - \vec{y}_1\|_2$', 'interpreter','latex')
ylabel('Analytical distance contributions')
h = legend('$\|\vec{y}_2 - \vec{y}_1\|_M^2$', '$\| $E$_M\|$',  '$ \|\vec{y}_2 - \vec{y}_1\|_M^2 + \|$E$_M \|$', 'location','southeast');
set(h,'Interpreter','latex');
set(h,'fontsize', 6);
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

%%
make_fig(3, 3.5);
loglog(dy2, Dis_avg, '.')
hold on
loglog(dy2, 0.5*dy2.^2, '-b')
for i=1:length(kernel_samplepoints)
    loglog(dy, kernel_samplepoints(i)*ones(size(dy)), 'linestyle', ':', 'color', 0.5*ones(1,3))
end
axis([1e-2 1e1 1e-6 1e2])
axis square
xlabel('$\| \vec{y}_2 - \vec{y}_1 \|_2$', 'interpreter','latex')
ylabel('$\| \vec{y}_2 - \vec{y}_1 \|_M^2$', 'interpreter','latex')
saveas(gcf, 'dist_dy_nonlinear.eps', 'epsc');

%%

dt_tmp = 10.^(-8:-2);
dt_tmp2 = logspace(-9, 0, 100);

make_fig(3, 3);
loglog(dt_tmp2, (1+4/epsilon+1/epsilon)*ones(size(dt_tmp2)), '-b')
hold on
% loglog(dt_tmp2, 8*sqrt(dt_tmp2)/epsilon^(3/2), '-r')
loglog(dt_tmp2, 15*dt_tmp2/epsilon^2, '-r')
% loglog(dt_tmp2, (1+4/epsilon+1/epsilon)+8*sqrt(dt_tmp2)/epsilon^(3/2), '-k')
loglog(dt_tmp2, (1+4/epsilon+1/epsilon)+15*dt_tmp2/epsilon^2, '-k')
for i=1:length(dt_burst_samplepoints)
    loglog(dt_burst_samplepoints(i)*ones(size(dt_tmp2)), 8*sqrt(dt_tmp2)/epsilon^(3/2), 'linestyle', ':', 'color', 0.5*ones(1,3))
end
axis([1e-9 1e-2 1e1 1e5])
axis square
xlabel('\delta t')
ylabel('Analytical covariance contributions')
h = legend('$\|C\|$', '$\| $E$_C\|$', '$\| C\| + \|$E$_C \|$', 'location','southeast');
set(h,'Interpreter','latex');
set(h,'fontsize', 6);
saveas(gcf, 'C_dt_analytical_nonlinear.eps', 'epsc');

return

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
make_fig(3.5, 3);
loglog(dt_tmp, mean(norm_c, 2),'.')
hold on
for i=1:length(dt_burst_samplepoints)
    loglog(dt_burst_samplepoints(i)*ones(size(dt_tmp2)), dt_tmp2/epsilon^2, 'linestyle', ':', 'color', 0.5*ones(1,3))
end
axis([1e-9 1e-2 1e1 1e5])
axis square
xlabel('\delta t')
ylabel('$\| \hat{C} \|$', 'interpreter','latex')
saveas(gcf, 'C_dt_nonlinear.eps', 'epsc');

%%

for j=1:length(dt_burst_samplepoints)
    
    dt_burst = dt_burst_samplepoints(j);
    
    [data_init, data_burst_init, t] = simulate_SDE(nsteps, dt, nsteps_burst, dt_burst, epsilon, data_init);
    
    data1 = f(data_init);
    data1_burst = zeros(size(data_burst_init));
    for i=1:nsteps
        data1_burst(:,:,i) = f(data_burst_init(:,:,i));
    end
    
    [c, inv_c, ranks] = covariances(data1_burst, dim, dt_burst);
    
    for j2=1:length(kernel_samplepoints)
        eps_NIV = kernel_samplepoints(j2);
        
        [V, D, eps, Dis] = NIV_return_dist(data1, inv_c, eps_NIV, 10, 0);
        if j==1 && j2==1
            V(:,2) = sign(corr(V(:,2), data_init(:,1)))*V(:,2);
        else
            V(:,2) = sign(corr(V(:,2), data_init(:,2)))*V(:,2);
        end
        
        make_fig(3, 3);
        scatter(data1(:,1),data1(:,2),50,V(:,2), '.')
        axis tight
        axis square
        xlabel('y_1')
        ylabel('y_2')
        h = colorbar;
        xlabel(h,'\phi_1');
        saveas(gcf, sprintf('data_nonlinear_NIV_dt%d_kernel%d.eps', j, j2), 'epsc');
    end
    
end




