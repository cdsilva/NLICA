clear all
close all

rng(123);
markersize = 500;


%% simulations


% rate of switching velocity
lambda = 400;
% speed of particles
s = sqrt(lambda);

% number of particles
N = 1000;
% maximum simulation time
tmax = 10;
tmin = 0;
% time step
dt = 1;
nsteps = floor(tmax / dt) + 1;

%initial probility of a particle moving to the right (vel=+1)
nsim = 10;
pinit = linspace(0.1, 0.9, nsim);

% histogram parameters
% nbins = 512;
nbins = 32;
x_hist = linspace(-s*tmax, s*tmax, nbins);
hist_all = zeros(nbins, nsteps*nsim);
hist_left = zeros(nbins, nsteps*nsim);
hist_right = zeros(nbins, nsteps*nsim);
all_time = zeros(1, nsteps*nsim);
all_p = zeros(1, nsteps*nsim);

%%
subplot_dim1 = floor(sqrt(nsim));
subplot_dim2 = ceil(nsim / subplot_dim1);

figure;
for j=1:nsim
    % simulate velocity jump process
    [pos, vel, time] = simulate_jumps(N, tmax, dt, lambda, s, pinit(j));
    %[pos, vel, time] = simulate_jumps_unifdist(N, tmax, dt, lambda, s, pinit(j), min_start, max_start);
    
    % histogram of all particle positions as a function of time
    for i=1:nsteps
        
        [hist_all(:,(j-1)*nsteps+i), ~] = hist(pos(:,i),x_hist);
        
        idx = (vel(:,i) == -1);
        [hist_left(:,(j-1)*nsteps+i), ~] = hist(pos(idx,i),x_hist);
        
        idx = (vel(:,i) == 1);
        [hist_right(:,(j-1)*nsteps+i), ~] = hist(pos(idx,i),x_hist);
        
    end
    
    % store times and probabilities
    all_time((j-1)*nsteps+1:j*nsteps) = time;
    all_p((j-1)*nsteps+1:j*nsteps) = pinit(j);
    
    % plot histograms
    %     subplot(subplot_dim1, subplot_dim2, j)
    %     plot(hist_all(:,(j-1)*nsteps+1:j*nsteps))
    %
end


%% dmaps on raw histograms

% use data after initial "relaxation" time
idx = (all_time > tmin);


%%

rho_diff = sum((hist_left-hist_right).^2);
rho_mean_left = x_hist * hist_left;
rho_mean_right = x_hist * hist_right;
rho_mean = x_hist * hist_all;
rho_mean_der = x_hist * (hist_all(:, 2:end) - hist_all(:, 1:end-1));
rho_mean_der = [0 rho_mean_der];

rho_mean_left_norm = x_hist * (hist_left./repmat(sum(hist_left), nbins, 1));
rho_mean_right_norm = x_hist * (hist_right./repmat(sum(hist_right), nbins, 1));

%% compute EMD between histograms


W2 = zeros(nsteps*nsim);

D_for_emd = pdist2(x_hist', x_hist');

for i1=1:nsteps*nsim
    for i2=1:i1-1
        
        %[emd_dist, ~]= emd_hat_gd_metric_mex(hist_all(:,i1), hist_all(:,i2), D_for_emd);
        emd_dist = sum(abs(cumsum(hist_all(:,i1)) - cumsum(hist_all(:,i2))));
        
        W2(i1, i2) = emd_dist^2;
        W2(i2, i1) = W2(i1, i2);
    end
end


%% dmaps on EMD

W2 = W2(idx,idx);

eps2 = median(W2(:));

[V2, D2] = dmaps(W2, eps2, 10);


if corr(V2(:,2), all_time(idx)') < 0
    V2(:,2) = -V2(:,2);
end
if corr(V2(:,3), all_p(idx)') < 0
    V2(:,3) = -V2(:,3);
end

%%
figure;
scatter(V2(:,2),V2(:,3),markersize,all_time(idx), '.')
xlabel('$\phi_1$', 'interpreter','latex')
ylabel('$\phi_2$', 'interpreter','latex')
h = colorbar;
set(get(h,'xlabel'),'String', 't');
%title(sprintf('EMD: colored by time, \\lambda= %d', lambda))
%print(sprintf('EMD_t_%d', lambda), '-r300','-djpeg')

figure;
scatter(V2(:,2),V2(:,3),markersize,all_p(idx), '.')
xlabel('$\phi_1$', 'interpreter','latex')
ylabel('$\phi_2$', 'interpreter','latex')
h = colorbar;
set(get(h,'xlabel'),'String', 'p');
%title(sprintf('EMD: colored by p, \\lambda= %d', lambda))
%print(sprintf('EMD_p_%d', lambda), '-r300','-djpeg')

corr(V2(:,2),all_time(idx)')

corr(V2(:,3),all_p(idx)')

return

%%

% figure;
% scatter(V2(:,2),V2(:,3),markersize,rho_diff(idx), '.')
% xlabel('$\phi_1$', 'interpreter','latex')
% ylabel('$\phi_2$', 'interpreter','latex')
% h = colorbar;
% set(get(h,'xlabel'),'String', 'FLUXES');
% 
% figure;
% scatter(V2(:,2),V2(:,3),markersize,rho_mean_left(idx), '.')
% xlabel('$\phi_1$', 'interpreter','latex')
% ylabel('$\phi_2$', 'interpreter','latex')
% h = colorbar;
% set(get(h,'xlabel'),'String', 'mean of left particles');
% 
% figure;
% scatter(V2(:,2),V2(:,3),markersize,rho_mean_right(idx), '.')
% xlabel('$\phi_1$', 'interpreter','latex')
% ylabel('$\phi_2$', 'interpreter','latex')
% h = colorbar;
% set(get(h,'xlabel'),'String', 'mean of right particles');

% figure;
% scatter(V2(:,2),V2(:,3),markersize,rho_mean(idx), '.')
% xlabel('$\phi_1$', 'interpreter','latex')
% ylabel('$\phi_2$', 'interpreter','latex')
% h = colorbar;
% set(get(h,'xlabel'),'String', 'mean of particles');
% 
% figure;
% scatter(V2(:,2),V2(:,3),markersize,rho_mean_der(idx), '.')
% xlabel('$\phi_1$', 'interpreter','latex')
% ylabel('$\phi_2$', 'interpreter','latex')
% h = colorbar;
% set(get(h,'xlabel'),'String', 'mean of derivative');

figure;
scatter(V2(:,2),V2(:,3),markersize,rho_mean_right(idx)-rho_mean_left(idx), '.')
xlabel('$\phi_1$', 'interpreter','latex')
ylabel('$\phi_2$', 'interpreter','latex')
h = colorbar;
set(get(h,'xlabel'),'String', 'E\rho^+ - E\rho^-');

figure;
scatter(V2(:,2),V2(:,3),markersize,rho_mean_right_norm(idx)-rho_mean_left_norm(idx), '.')
xlabel('$\phi_1$', 'interpreter','latex')
ylabel('$\phi_2$', 'interpreter','latex')
h = colorbar;
set(get(h,'xlabel'),'String', 'E\rho^+_{norm} - E\rho^-_{norm}');