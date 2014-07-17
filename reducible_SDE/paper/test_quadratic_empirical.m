clear all
close all


%% define parameters

dim = 2;

epsilon = 1e-3;


dt = 1e-4;
dt_burst = 1e-10;

nsteps = 3000;
nsteps_burst = 200;


%%

dy_axis_lim = [1e-2 1e1 1e-10 1e5];

make_fig;
deltas = logspace(-2, 1, 100);
loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2)
hold on
loglog(deltas, (38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2, '-g')
loglog(deltas, 10*deltas.^4, '-r')
loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2+(38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4, '-k')
% loglog(sqrt(eps_dmaps/(0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2)))), eps_dmaps, 'ok')
% legend('linear approximation distance','error from covariance estimation','error from Taylor expansion','location','best')
xlabel('$\| Y_2 - Y_1 \|$','interpreter','latex')
ylabel('$\| X_2 - X_1 \|^2_{est}$','interpreter','latex')
axis(dy_axis_lim)
print(gcf, '-depsc', 'errors_function_dy');

% make_fig;
% loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2+(38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4, '-k')
% hold on
% loglog(deltas, 2*deltas.^2, '-r')
% xlabel('$\| Y_2 - Y_1 \|$','interpreter','latex')
% ylabel('$\| X_2 - X_1 \|^2_{est}$','interpreter','latex')
% axis(dy_axis_lim)
% print(gcf, '-depsc', 'totaldist_function_dy');

%% simulate SDE

rng(123);

[data_init, data_burst_init, data1, data1_burst, t] = simulate_quadratic(nsteps, dt, nsteps_burst, dt_burst, epsilon);

% idx = find(data_init(:, 2).^2 < epsilon/2);
% data_init = data_init(idx, :);
% t = t(idx, :);
% nsteps = length(idx);

make_fig;
scatter(data_init(:,1),data_init(:,2),50,t,'.')
xlabel('X(1)')
ylabel('X(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal

make_fig;
scatter(data1(:,1),data1(:,2),50,t,'.')
xlabel('Y(1)')
ylabel('Y(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal

%%

[inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);
inv_c1 = inv_c1 * dt_burst;

%%

[~, ~, ~, Dis_data1] = NIV_return_dist(data1, inv_c1, 0, 10, 0);

Dis_Y = squareform(pdist(data1));

idx = floor(linspace(1,numel(Dis_Y), 1000));
make_fig;
loglog(Dis_Y(idx), Dis_data1(idx), '.')
hold on
% loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2+(38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4, '-r')
loglog(deltas(1:end-1), 0.3*deltas(1:end-1).^2, '-r')
xlabel('$\| Y_2 - Y_1 \|$', 'interpreter','latex')
ylabel('$\| X_2 - X_1 \|^2_{est}$', 'interpreter','latex')
axis(dy_axis_lim)
print(gcf, '-depsc', 'empirical_totaldist_function_dy');

%%
deltas = logspace(-2, 0.8, 20);
curve_means = zeros(size(deltas));
curve_means = curve_means(1:end-1);
[bincounts, ind] = histc(Dis_Y(:), deltas);
for i=1:length(deltas)-1
    curve_means(i) = mean(Dis_data1(ind == i));
end

make_fig;
loglog(deltas(1:end-1), curve_means, '.b', 'markersize', 5)
hold on
loglog(deltas(1:end-1), deltas(1:end-1).^2, '-r')
hold on
% loglog(deltas, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2+(38 * dt_burst/(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4, '-k')
xlabel('$\| Y_2 - Y_1 \|$', 'interpreter','latex')
ylabel('$\| X_2 - X_1 \|^2_{est}$', 'interpreter','latex')
axis(dy_axis_lim)
% print(gcf, '-depsc', 'empirical_totaldist_function_dy');

%%

dt_axis_lim = [1e-10 1e0 1e-8 1e-2];

make_fig;
deltas = 0.01;
dt_burst = logspace(-10, -2, 100);
loglog(dt_burst, 0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2*ones(size(dt_burst)))
hold on
loglog(dt_burst, (38 * dt_burst./(epsilon^2 + 2 * dt_burst))*deltas.^2, '-g')
loglog(dt_burst, 10*deltas.^4*ones(size(dt_burst)), '-r')
loglog(dt_burst,  0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2*ones(size(dt_burst))+(38 * dt_burst./(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4*ones(size(dt_burst)), '-k')
axis(dt_axis_lim)
xlabel('$\delta t$','interpreter','latex')
ylabel('$\| X_2 - X_1 \|^2_{est}$','interpreter','latex')
print(gcf, '-depsc', 'errors_function_dt');

make_fig;
loglog(dt_burst,  0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2*ones(size(dt_burst))+(38 * dt_burst./(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4*ones(size(dt_burst)), '-k')
axis(dt_axis_lim)
xlabel('$\delta t$','interpreter','latex')
ylabel('$\| X_2 - X_1 \|^2_{est}$','interpreter','latex')
print(gcf, '-depsc', 'totaldist_function_dt');

%%

dt_burst = 10.^(-10:-2);
mean_Dis = zeros(size(dt_burst));

nsteps = 2;

% data_init = [0 0;
%     0 deltas*sqrt(0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2)))];
data_init = [0 0;
    deltas*sqrt(0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))) 0];

for j=1:length(dt_burst)
    
    [~, data_burst_init, data1, data1_burst, ~] = simulate_quadratic(0, 0, nsteps_burst, dt_burst(j), epsilon, data_init);
    
    [inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);
    inv_c1 = inv_c1 * dt_burst(j);

    [~, ~, ~, Dis_data1] = NIV_return_dist(data1, inv_c1, 0, 1, 0);

    mean_Dis(j) = Dis_data1(1,2);
end

make_fig;
loglog(dt_burst, mean_Dis, '.')
hold on
% loglog(dt_burst,  0.5*(2 + epsilon + sqrt(9 + 2 * epsilon + epsilon^2))*deltas.^2*ones(size(dt_burst))+(38 * dt_burst./(epsilon^2 + 2 * dt_burst))*deltas.^2+10*deltas.^4*ones(size(dt_burst)), '-k')
xlabel('$\delta t$','interpreter','latex')
ylabel('$\| X_2 - X_1 \|^2_{est}$','interpreter','latex')
axis(dt_axis_lim)
print(gcf, '-depsc', 'empirical_totaldist_function_dt');

%%

eps_dmaps = 0.01;
dt_burst = 1e-9;
dt = 1e-4;


nsteps = 3000;
nsteps_burst = 200;


rng(123);
[data_init, data_burst_init, data1, data1_burst, t] = simulate_quadratic(nsteps, dt, nsteps_burst, dt_burst, epsilon);


make_fig;
scatter(data_init(:,1),data_init(:,2),50,t,'.')
xlabel('X(1)')
ylabel('X(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal
print(gcf, '-depsc', 'original_data');

make_fig;
scatter(data1(:,1),data1(:,2),50,t,'.')
xlabel('Y(1)')
ylabel('Y(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal
print(gcf, '-depsc', 'transformed_data');


%% DMAPS

W = squareform(pdist(data_init)).^2;
[V_init_dmaps, D_init_dmaps] = dmaps(W, eps_dmaps, 10);

make_fig;
scatter(data_init(:,1),data_init(:,2),50,V_init_dmaps(:,2),'.')
xlabel('x')
ylabel('y')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal


%%

[inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);
inv_c1 = inv_c1 * dt_burst;

%%

eps_niv = eps_dmaps;

[V1_niv, D1_niv, ~, Dis_data1] = NIV_return_dist(data1, inv_c1, eps_niv, 10, 0);


make_fig;
scatter(data_init(:,1),data_init(:,2),50,V1_niv(:,2),'.')
xlabel('x')
ylabel('y')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal

make_fig;
plot(V1_niv(:,2), V_init_dmaps(:,2),'.')
xlabel('NIV from Y data')
ylabel('DMAPS from X data')
axis equal
print(gcf, '-depsc', 'NIV_correlation');

