clear all
close all


%% define parameters

dim = 2;

epsilon = 1e-3;


%%
    
eps_dmaps = 0.01;
dt_burst = 1e-9;
dt = 1e-4;


nsteps = 3000;
nsteps_burst = 200;


rng(123);
[data_init, data_burst_init, data1, data1_burst, t] = simulate_quadratic(nsteps, dt, nsteps_burst, dt_burst, epsilon);


make_fig;
scatter(data_init(:,1),data_init(:,2),50,t*dt,'.')
xlabel('X(1)')
ylabel('X(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal
print(gcf, '-depsc', 'data_init');


make_fig;
scatter(data1(:,1),data1(:,2),50,t*dt,'.')
xlabel('Y(1)')
ylabel('Y(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', 't');
axis equal
print(gcf, '-depsc', 'data_transformed');


%% DMAPS

W = squareform(pdist(data_init)).^2;
[V_init_dmaps, D_init_dmaps] = dmaps(W, eps_dmaps, 10);


%%

[inv_c1, new_data1, ranks1] = covariances2(data1_burst, dim);
inv_c1 = inv_c1 * dt_burst;

%%

eps_niv = eps_dmaps;

[V1_niv, D1_niv, ~, Dis_data1] = NIV_return_dist(data1, inv_c1, eps_niv, 10, 0);


make_fig;
scatter(data_init(:,1),data_init(:,2),50,V1_niv(:,2),'.')
xlabel('X(1)')
ylabel('X(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal

make_fig;
scatter(data1(:,1),data1(:,2),50,V1_niv(:,2),'.')
xlabel('Y(1)')
ylabel('Y(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal
print(gcf, '-depsc', 'data_transformed_colored_NIV');


make_fig;
plot(V1_niv(:,2), V_init_dmaps(:,2),'.')
xlabel('NIV from Y data')
ylabel('DMAPS from X data')
axis equal


make_fig;
plot(data_init(:,1), V1_niv(:,2),'.')
xlabel('X(1)')
ylabel('NIV from Y data')


%%

W = squareform(pdist(data1)).^2;
eps1 = max(median(W));
[V1_dmaps, D1_dmaps] = dmaps(W, eps1, 10);

make_fig;
scatter(data1(:,1),data1(:,2),50,V1_dmaps(:,2),'.')
xlabel('Y(1)')
ylabel('Y(2)')
h = colorbar('peer',gca);
set(get(h,'xlabel'),'String', '\phi_1');
axis equal
print(gcf, '-depsc', 'data_transformed_colored_DMAPS');

make_fig;
plot(data_init(:,1), V1_dmaps(:,2),'.')
xlabel('X(1)')
ylabel('NIV from Y data')


