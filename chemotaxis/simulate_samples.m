function W2 = simulate_samples(lambda, s)


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

% histogram parameters
nbins = 32;
x_hist = linspace(-s*tmax, s*tmax, nbins);
hist_all = zeros(nbins, nsteps*nsim);
hist_left = zeros(nbins, nsteps*nsim);
hist_right = zeros(nbins, nsteps*nsim);
all_time = zeros(1, nsteps*nsim);
all_p = zeros(1, nsteps*nsim);

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
        
    end
    
    % store times and probabilities
    all_time((j-1)*nsteps+1:j*nsteps) = time;
    all_p((j-1)*nsteps+1:j*nsteps) = pinit(j);
    
end

% use data after initial "relaxation" time
idx = (all_time > 1);

%% compute EMD

W2 = zeros(nsteps*nsim);

for i1=1:nsteps*nsim
    for i2=1:i1-1
        
        %[emd_dist, ~]= emd_hat_gd_metric_mex(hist_all(:,i1), hist_all(:,i2), D_for_emd);
        emd_dist = sum(abs(cumsum(hist_all(:,i1)) - cumsum(hist_all(:,i2))));
        
        W2(i1, i2) = emd_dist^2;
        W2(i2, i1) = W2(i1, i2);
    end
end

W2 = W2(idx,idx);

