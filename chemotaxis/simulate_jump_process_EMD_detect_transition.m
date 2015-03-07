clear all
close all

rng(123);
markersize = 500;


%% simulations

nsim2 = 10;

lambda_all = linspace(1, 1000, nsim2);
s_all = sqrt(lambda_all);

% number of particles
N = 1000;
% time step
dt = 0.1;
% maximum simulation time
tmax = 10;
nsteps = floor(tmax / dt) + 1;

%initial probility of a particle moving to the right (vel=+1)
nsim = 10;
pinit = linspace(0.1, 0.9, nsim);

min_start = -10;
max_start = 10;

% histogram parameters
%nbins = 128;
nbins = 32;
hist_all = zeros(nbins, nsteps*nsim);
hist_left = zeros(nbins, nsteps*nsim);
hist_right = zeros(nbins, nsteps*nsim);
all_time = zeros(1, nsteps*nsim);
all_p = zeros(1, nsteps*nsim);

eval1_vec = [];
eval2_vec = [];
t_vec = [];
lam_vec = [];

% f1 = figure;
% f2 = figure;
% fig_idx = 1;

niter = 3;

for k=1:niter
    for sim_num = 1:length(lambda_all)
        
        sim_num
        
        lambda = lambda_all(sim_num);
        s = s_all(sim_num);
        
        x_hist = linspace(min_start-s*tmax, max_start+s*tmax, nbins);
        
        
        %%
        for j=1:nsim
            % simulate velocity jump process
            [pos, vel, time] = simulate_jumps(N, tmax, dt, lambda, s, pinit(j));
            
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
        
        
        
        %% compute EMD between histograms
        
        W = zeros(nsteps*nsim);
        
        D_for_emd = pdist2(x_hist', x_hist');
        
        for i1=1:nsteps*nsim
            for i2=1:i1-1
                
                emd_dist = sum(abs(cumsum(hist_all(:,i1)) - cumsum(hist_all(:,i2))));
                
                W(i1, i2) = emd_dist^2;
                W(i2, i1) = W(i1, i2);
            end
        end
        
        %% dmaps on EMD
        
        
        for i1=20*dt:10*dt:tmax
            idx = find(all_time > i1-10*dt & all_time <= i1);
            W2 = W(idx,idx);
            eps2 = median(W2(:));
            [V2, D2] = dmaps(W2, eps2, 7);
            
            %
            eps_med_scale = 3;
            res = compute_residuals_DMAPS(V2, eps_med_scale);
            
            
            %         i2 = find(res > 0.3);
            [~, i2] = sort(res, 'descend');
            i2 = sort(i2(1:2));
            
            eval1_vec = [eval1_vec; D2(i2(1), i2(1))];
            eval2_vec = [eval2_vec; D2(i2(2), i2(2))];
            t_vec = [t_vec; i1];
            lam_vec = [lam_vec; lambda];
            
%                     figure(f1);
%                     subplot(nsim2, 10, fig_idx)
%                     scatter(V2(:, i2(1)), V2(:, i2(2)),50, all_p(idx), '.')
%                     title(sprintf('t=%d, \\lambda=%4.0f', i1, lambda))
%             %         axis off
%             
%                     figure(f2);
%                     subplot(nsim2, 10, fig_idx)
%                     scatter(V2(:, i2(1)), V2(:, i2(2)),50, all_time(idx), '.')
%                     title(sprintf('t=%d, \\lambda=%4.0f', i1, lambda))
%             %         axis off
%             
%                     fig_idx = fig_idx + 1;
            
        end
        
    end
end


%%

X = unique(lam_vec);
Y = unique(t_vec);
C = zeros(length(X), length(Y));
C2 = zeros(length(X), length(Y));

for i=1:length(eval1_vec)
    i1 = find(X == lam_vec(i));
    i2 = find(Y == t_vec(i));
    C(i1, i2) = C(i1, i2) + sqrt(log(eval1_vec(i))./log(eval2_vec(i)));
    C2(i1, i2) = C2(i1, i2) + max(lam_vec(i)*t_vec(i), 1/(lam_vec(i)*t_vec(i)));
end

% C  = C / niter;
% C2  = C2 / niter;


figure;
imagesc(X, Y, C')
set(gca, 'xdir','normal')
set(gca, 'ydir','normal')
colorbar


figure;
imagesc(X, Y, (C2'))
set(gca, 'xdir','normal')
set(gca, 'ydir','normal')
colorbar