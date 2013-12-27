clear all
close all

%% simulations

% rate of switching velocity
lambda = 1;
% speed of particles
s = 10;

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
nbins = 128;
x_hist = linspace(min_start-s*tmax, max_start+s*tmax, nbins);
hist_all = zeros(nbins, nsteps*nsim);
hist_left = zeros(nbins, nsteps*nsim);
hist_right = zeros(nbins, nsteps*nsim);
all_time = zeros(1, nsteps*nsim);
all_p = zeros(1, nsteps*nsim);

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
    subplot(subplot_dim1, subplot_dim2, j)
    plot(hist_all(:,(j-1)*nsteps+1:j*nsteps))
    
end

%% dmaps on raw histograms

% use data after initial "relaxation" time
idx = (all_time > 1);

W = squareform(pdist(hist_all(:,idx)')).^2;
eps = median(W(:));

[V, D] = dmaps(W, eps, 10);

figure;
scatter(V(:,2),V(:,3),50,all_time(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('histograms: colored by time')

figure;
scatter(V(:,2),V(:,3),50,all_p(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('histograms: colored by p')

%% compute scattering coefficients for histograms

%addpath 'C:\Users\cdsilva\Documents\MATLAB\scatnet-0.2';
%addpath_scatnet

filt_opt = struct();
filt_opt.filter_type = 'morlet_1d';
%filt_opt.Q = 1;
%filt_opt.J = T_to_J(1024, filt_opt.Q); % 1024 is the analysis window size.

sc_opt = struct();
sc_opt.antialiasing = 1; % by default 1 -> 2^1

cascade = wavelet_factory_1d(nbins, filt_opt, sc_opt); 

% store scattering coefficients
S_all = [];

for i=1:nsteps*nsim
    % compute scattering coefficients for each histogram
    [S, U] = scat(hist_all(:,i), cascade);

    %Stilde = average_scat(S,10*2^10);
    %U{1}.meta.resolution = 0;
    %Utilde = average_scat(U,10*2^10);

    %S = log_scat(renorm_scat(S));

    [S_table, meta] = format_scat(S);
    %S_table = reshape(S_table, [size(S_table,1), size(S_table,3)]);

    % store scattering coefficients
    S_all = [S_all S_table(:)];
end

%% dmaps on scattering coefficients
W2 = squareform(pdist(S_all(:,idx)')).^2;
eps2 = median(W2(:));

[V2, D2] = dmaps(W2, eps2, 10);

figure;
scatter(V2(:,2),V2(:,3),50,all_time(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('scattering transform: colored by time')

figure;
scatter(V2(:,2),V2(:,3),50,all_p(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('scattering transform: colored by p')

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

figure;
scatter(V2(:,2),V2(:,3),50,all_time(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('EMD: colored by time')

figure;
scatter(V2(:,2),V2(:,3),50,all_p(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('EMD: colored by p')


