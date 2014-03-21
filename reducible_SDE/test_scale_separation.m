clear all
close all

a = 5;
b = 1;

f = @(x) a;
g = @(x) -b*x;

eps = 0.0001;

drift = @(t, x) [f(x(1)); g(x(2))/eps];
diffn = @(t, x) [1 0; 0 1/sqrt(eps)];

x0 = [0; 0];
tmax = 0.1;

dt = 0.00001;
nsteps_per_step = 10;

SDE = sde(drift, diffn, 'StartState', x0);
nPeriods = ceil(tmax / dt);
[x, t, Z] = SDE.simulate(nPeriods, 'DeltaTime', dt, 'NSTEPS', nsteps_per_step);

figure; scatter(x(1:5:end,1),x(1:5:end,2), 50, t(1:5:end),'.')
xlabel('x_1')
ylabel('x_2')
title('Colored by time')
print('raw_data.jpg', '-djpeg', '-r300')

%%
t_window = 0.001;
knn = round((t_window/2) / dt);
D = 2;
stride = 4;
[inv_c, new_x, ranks] = covariances(x, knn, D, stride);
%inv_c = inv_c * dt1;

%% NIV
if size(new_x,1) > 4000
    disp('Too much data; NIV will fail')
    return
end

neigs = 10;
[V, D, ~] = NIV(new_x, inv_c, 0, neigs, 0);
V = decouple_dmaps(V);


%%
figure;
scatter(new_x(:,1),new_x(:,2),200,V(:,2),'.')
xlabel('x_1')
ylabel('x_2')
title('NIV: Colored by \phi_2')
print('NIV_data.jpg', '-djpeg', '-r300')

figure;
plot(new_x(:,1),V(:,2),'.')
xlabel('x_1')
ylabel('\phi_2 from NIV')
print('NIV_corr.jpg', '-djpeg', '-r300')

%%

W = squareform(pdist(new_x)).^2;
eps2 = median(W(:));
[V2, D2] = dmaps(W, eps2, neigs);

figure;
scatter(new_x(:,1),new_x(:,2),200,V2(:,2),'.')
xlabel('x_1')
ylabel('x_2')
title('DMAPS: Colored by \phi_2')
print('DMAPS_data.jpg', '-djpeg', '-r300')

figure;
plot(new_x(:,1),V2(:,2),'.')
xlabel('x_1')
ylabel('\phi_2 from DMAPS')
print('DMAPS_corr.jpg', '-djpeg', '-r300')

%%

figure;
plot(V2(:,2),V2(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')
title('DMAPS')
print('DMAPS_embed.jpg', '-djpeg', '-r300')

figure;
plot(V(:,2),V(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')
title('NIV')
print('NIV_embed.jpg', '-djpeg', '-r300')

