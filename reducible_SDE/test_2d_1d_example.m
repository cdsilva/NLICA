clear all
close all

rng(123);

b = 1;
a = @(x) -100*b*x^99;

drift = @(t, x) [a(x(1)); a(x(2))];
diffn = @(t, x) [1 0; 0 1];

theta0 = [0; 0];
tmax = 10;

dim = 2;
subsample_step = 80;

dt = 0.0001;
nsteps_per_step = 10;

SDE = sde(drift, diffn, 'StartState', theta0);
nPeriods = ceil(tmax / dt);
[theta, t, Z] = SDE.simulate(nPeriods, 'DeltaTime', 4*dt, 'NSTEPS', nsteps_per_step);

theta = (theta+1)/2;

figure; scatter(theta(:,1),theta(:,2), 50, t,'.')
xlabel('x_1')
ylabel('x_2')
title('Colored by time')

%%
window_time = 0.1;
window_steps = round(window_time/2 / dt);

start_idx = 60;
end_idx = start_idx+window_steps;

figure; 
hist(theta(start_idx+1:end_idx,1)-theta(start_idx:end_idx-1,1),20)

figure; plot(theta(:,1),theta(:,2),'.')
hold on
plot(theta(start_idx:end_idx,1),theta(start_idx:end_idx,2),'.k', 'markersize', 10)



%%
%f = @(x) [(x(:,1)+1)/2+((x(:,2)+1)/2).^3 (x(:,2)+1)/2-((x(:,1)+1)/2).^3 zeros(size(x,1), 1)];

r = 0.1;
f = @(x) [r*(x(:,1)+x(:,2).^3)./sqrt((x(:,1)+x(:,2).^3).^2+((x(:,2)-x(:,1).^3).^2)+1) r*(x(:,2)-x(:,1).^3)./sqrt((x(:,1)+x(:,2).^3).^2+((x(:,2)-x(:,1).^3).^2)+1) r./sqrt((x(:,1)+x(:,2).^3).^2+((x(:,2)-x(:,1).^3).^2)+1)];


J = @(x) r*[(1-2*x(:,1).^6+x(:,1).^3.*x(:,2)+x(:,2).^2-3*x(:,1).^5.*x(:,2).^3+3*x(:,1).^2.*x(:,2).^4)./(1+(x(:,1).^3-x(:,2)).^2+(x(:,1)+x(:,2).^3).^2).^(3/2) (x(:,1).^4-x(:,1).*x(:,2)+3*x(:,2).^2+3*x(:,1).^6.*x(:,2).^2-5*x(:,1).^3.*x(:,2).^3+2*x(:,2).^4)./(1+(x(:,1).^3-x(:,2)).^2+(x(:,1)+x(:,2).^3).^2).^(3/2) ...
    (-2*x(:,1).^4-x(:,1).*x(:,2)-5*x(:,1).^3.*x(:,2).^3-x(:,2).^4-3*x(:,1).^2.*(1+x(:,2).^6))./(1+(x(:,1).^3-x(:,2)).^2+(x(:,1)+x(:,2).^3).^2).^(3/2) (1+x(:,1).^2+3*x(:,1).^4.*x(:,2).^2-x(:,1).*x(:,2).^3+3*x(:,1).^3.*x(:,2).^5-2*x(:,2).^6)./(1+(x(:,1).^3-x(:,2)).^2+(x(:,1)+x(:,2).^3).^2).^(3/2) ...
    (-x(:,1)-3*x(:,1).^5+3*x(:,1).^2.*x(:,2)-x(:,2).^3)./(1+(x(:,1).^3-x(:,2)).^2+(x(:,1)+x(:,2).^3).^2).^(3/2) (x(:,1).^3-x(:,2)-3*x(:,1).*x(:,2).^2-3*x(:,2).^5)./(1+(x(:,1).^3-x(:,2)).^2+(x(:,1)+x(:,2).^3).^2).^(3/2)];

x = f(theta);
figure;
scatter3(x(1:subsample_step:end,1),x(1:subsample_step:end,2),x(1:subsample_step:end,3),50, t(1:subsample_step:end), '.')

%%

drift_noise = @(t, x) [a(x(1)); a(x(2)); a(x(3))];
diffn_noise = @(t, x) eye(3);

theta0 = [0; 0; 0];

%tscale_noise = 1000;
tscale_noise = 10;

dt2 = dt * tscale_noise;

nsteps_per_step = 10;

SDE = sde(drift_noise, diffn_noise, 'StartState', theta0);
[xi, ~, ~] = SDE.simulate(nPeriods, 'DeltaTime', dt2, 'NSTEPS', nsteps_per_step);

%%

[inv_c_true, ~, ~] = covariances_diff(x, window_steps, dim, subsample_step);

covar_diff = [];

figure;
lscale_noise = [0, 0.01, 0.1,1];
for i=1:4;
x_noise = x + xi * lscale_noise(i);

% figure;
% scatter3(x_noise(:,1),x_noise(:,2), x_noise(:,3),50, t, '.')
% 
% figure;
% scatter3(xi(:,1),xi(:,2), xi(:,3),50, t, '.')
% 
% figure; 
% hist(x(2:end,1)-x(1:end-1,1),100)
% 
% figure; 
% hist(x_noise(2:end,1)-x_noise(1:end-1,1),100)
% 


[inv_c, subsample_x, ranks] = covariances_diff(x_noise, window_steps, dim, subsample_step);
covar_diff = [covar_diff; norm(inv_c(:)-inv_c_true(:))];

if size(subsample_x,1) > 4000
    disp('Too much data; NIV will fail')
    return
end
neigs = 10;
eps_NIV = 0;
[V, D, eps_NIV] = NIV(subsample_x, inv_c, eps_NIV, neigs, 0);
V = decouple_dmaps(V);

subplot(2,2,i)
plot(V(:,2),V(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')
title(sprintf('size of noise =  %2.2f', lscale_noise(i)))
% figure;
% plot(theta(1:subsample_step:end,1),V(:,2),'.')
% 
% figure;
% plot(theta(1:subsample_step:end,2),V(:,3),'.')
end

figure; bar(covar_diff)
title('ERROR')

%%

W = squareform(pdist(subsample_x)).^2;
eps_DMAPS = median(W(:));

[V2, D2] = dmaps(W, eps_DMAPS, neigs);
V2 = decouple_dmaps(V2);

figure;
plot(V2(:,2),V2(:,3),'.')

figure;
plot(theta(1:subsample_step:end,1),V2(:,2),'.')

figure;
plot(theta(1:subsample_step:end,2),V2(:,3),'.')






