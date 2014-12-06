clear all
close all

rng(321);

markersize = 500;

%% simulations

% % rate of switching velocity
% lambda = 0.01;
% % speed of particles
% s = 0.1;

% % rate of switching velocity
% lambda = 1;
% % speed of particles
% s = 1;

% % rate of switching velocity
% lambda = 400;
% % speed of particles
% s = 20;

% % rate of switching velocity
% lambda = 2500;
% % speed of particles
% s = 50;

% % rate of switching velocity
% lambda = 4900;
% % speed of particles
% s = 70;


% number of particles
N = 1000;
% time step
dt = 1;
% maximum simulation time
tmax = 10;
nsteps = floor(tmax / dt) + 1;

%initial probility of a particle moving to the right (vel=+1)
nsim = 10;
pinit = linspace(0.1, 0.9, nsim);

min_start = -10;
max_start = 10;

% histogram parameters
nbins = 32;
hist_all = zeros(nbins, nsteps*nsim);
hist_left = zeros(nbins, nsteps*nsim);
hist_right = zeros(nbins, nsteps*nsim);
all_time = zeros(1, nsteps*nsim);
all_p = zeros(1, nsteps*nsim);
mean_pos = zeros(1, nsteps*nsim);
mean_vel = zeros(1, nsteps*nsim);

figure;

for lambda = [1 400 2500 4900]
    s = sqrt(lambda);
    
    x_hist = linspace(min_start-s*tmax, max_start+s*tmax, nbins);

%%

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
        
        mean_pos((j-1)*nsteps+i) = mean(pos(:,i));
        mean_vel((j-1)*nsteps+i) = mean(vel(:,i));
    end
    
    % store times and probabilities
    all_time((j-1)*nsteps+1:j*nsteps) = time;
    all_p((j-1)*nsteps+1:j*nsteps) = pinit(j);
    
end


% use data after initial "relaxation" time
idx = (all_time > 1);


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

if corr(V2(:,2), all_p(idx)') < 0
    V2(:,2) = -V2(:,2);
end
if corr(V2(:,3), all_time(idx)') < 0
    V2(:,3) = -V2(:,3);
end

%%
% figure;
% scatter(V2(:,2),V2(:,3),markersize,all_time(idx), '.')
% xlabel('$\phi_1$', 'interpreter','latex')
% ylabel('$\phi_2$', 'interpreter','latex')
% h = colorbar;
% set(get(h,'xlabel'),'String', 't');
% 
% figure;
% scatter(V2(:,2),V2(:,3),markersize,all_p(idx), '.')
% xlabel('$\phi_1$', 'interpreter','latex')
% ylabel('$\phi_2$', 'interpreter','latex')
% h = colorbar;
% set(get(h,'xlabel'),'String', 'p');

%%


%%
potential_cv = [all_p(idx)' all_time(idx)' mean_pos(idx)' mean_vel(idx)'];
[B, FitInfo] = lasso(potential_cv, V2(:,2));


figure; 
plot(B')
legend('p','t','pos','vel')

end

