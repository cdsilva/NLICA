clear all
close all

rng(321);

%% simulations

% rate of switching velocity
lambda = 100;
% speed of particles
s = 10;

% % rate of switching velocity
% lambda = 400;
% % speed of particles
% s = 20;

% number of particles
N = 1000;


% % time step
% dt = 2;
% % maximum simulation time
% tmax = 20;

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
x_hist = linspace(min_start-s*tmax, max_start+s*tmax, nbins);
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
idx = (all_time > 1);

% W = squareform(pdist(hist_all(:,idx)')).^2;
% eps = median(W(:));
% 
% [V, D] = dmaps(W, eps, 10);
% 
% figure;
% scatter(V(:,2),V(:,3),200,all_time(idx), '.')
% xlabel('$\phi_2$', 'interpreter','latex', 'fontsize', 20)
% ylabel('$\phi_3$', 'interpreter','latex', 'fontsize', 20)
% %title('histograms: colored by time')
% h = colorbar;
% set(get(h,'xlabel'),'String', 't', 'fontsize', 20);
% 
% figure;
% scatter(V(:,2),V(:,3),200,all_p(idx), '.')
% xlabel('$\phi_2$', 'interpreter','latex', 'fontsize', 20)
% ylabel('$\phi_3$', 'interpreter','latex', 'fontsize', 20)
% %title('histograms: colored by p')
% h = colorbar;
% set(get(h,'xlabel'),'String', 'p', 'fontsize', 20);

%% compute EMD between histograms

addpath('FastEMD');

W2 = zeros(nsteps*nsim);

D_for_emd = pdist2(x_hist', x_hist');

for i1=1:nsteps*nsim
    for i2=1:i1-1
        
        [emd_dist, ~]= emd_hat_gd_metric_mex(hist_all(:,i1), hist_all(:,i2), D_for_emd);
        
        
        W2(i1, i2) = emd_dist^2;
        W2(i2, i1) = W2(i1, i2);
    end
end

rmpath('FastEMD');

%% dmaps on EMD

W2 = W2(idx,idx);

eps2 = median(W2(:));

[V2, D2] = dmaps(W2, eps2, 10);

% if lambda <= 1
%     if corr(V2(:,2), all_p(idx)') < 0
%         V2(:,2) = -V2(:,2);
%     end
%     if corr(V2(:,3), all_time(idx)') < 0
%         V2(:,3) = -V2(:,3);
%     end
% end
% if lambda > 1
%     if corr(V2(:,2), all_time(idx)') < 0
%         V2(:,2) = -V2(:,2);
%     end
%     if corr(V2(:,3), all_p(idx)') < 0
%         V2(:,3) = -V2(:,3);
%     end
% end


%%
figure;
scatter(V2(:,2),V2(:,3),200,all_time(idx), '.')
xlabel('$\phi_2$', 'interpreter','latex')
ylabel('$\phi_3$', 'interpreter','latex')
h = colorbar;
set(get(h,'xlabel'),'String', 't');
%title(sprintf('EMD: colored by time, \\lambda= %d', lambda))
print(sprintf('EMD2_t_%d', tmax), '-r300','-djpeg')

figure;
scatter(V2(:,2),V2(:,3),200,all_p(idx), '.')
xlabel('$\phi_2$', 'interpreter','latex')
ylabel('$\phi_3$', 'interpreter','latex')
h = colorbar;
set(get(h,'xlabel'),'String', 'p');
%title(sprintf('EMD: colored by p, \\lambda= %d', lambda))
print(sprintf('EMD2_p_%d', tmax), '-r300','-djpeg')



