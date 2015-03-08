clear all
close all

rng(123);


%% simulations

nlambda = 10;

lambda_all = linspace(1, 1000, nlambda);
s_all = sqrt(lambda_all);

% number of particles
N = 1000;
% time step
dt = 1;
% maximum simulation time
tmax = 100;
nsteps = floor(tmax / dt) + 1;

t_range_all = 10:10:tmax;
ntime = length(t_range_all);

%initial probility of a particle moving to the right (vel=+1)
nsim = 10;
pinit = linspace(0.1, 0.9, nsim);

% histogram parameters
nbins = 256;
% nbins = 32;
hist_all = zeros(nbins, nsteps*nsim);
hist_left = zeros(nbins, nsteps*nsim);
hist_right = zeros(nbins, nsteps*nsim);
all_time = zeros(1, nsteps*nsim);
all_p = zeros(1, nsteps*nsim);

niter = 1;

eval1 = zeros(nlambda, ntime, niter);
eval2 = zeros(nlambda, ntime, niter);
t_dom = zeros(nlambda, ntime, niter);

f1 = figure;
f2 = figure;
fig_idx = 1;

%%

for k=1:niter
    for sim_num = 1:nlambda
        
        sim_num
        
        lambda = lambda_all(sim_num);
        s = s_all(sim_num);
        
        x_hist = linspace(-s*tmax, s*tmax, nbins);
        
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
        
        
        
        % compute EMD between histograms
        
        W = zeros(nsteps*nsim);
        
        for i1=1:nsteps*nsim
            for i2=1:i1-1
                
                emd_dist = sum(abs(cumsum(hist_all(:,i1)) - cumsum(hist_all(:,i2))));
                
                W(i1, i2) = emd_dist^2;
                W(i2, i1) = W(i1, i2);
            end
        end
        
        % dmaps on EMD
        
        for i1=1:ntime
            
            %% TODO check-- min time be 1??
            idx = find(all_time > t_range_all(i1)-9 & all_time <= t_range_all(i1));
            W2 = W(idx,idx);
            eps2 = median(W2(:));
            [V2, D2] = dmaps(W2, eps2, 10);
            
            %
            eps_med_scale = 3;
            res = compute_residuals_DMAPS(V2, eps_med_scale);
            
            %         i2 = find(res > 0.3);
            [~, i2] = sort(res, 'descend');
            i2 = sort(i2(1:2));
            %         i2 = [2 3];
            
            eval1(sim_num, i1, k) = D2(i2(1), i2(1));
            eval2(sim_num, i1, k) = D2(i2(2), i2(2));
            
            if abs(corr(all_time(idx)', V2(:, i2(1)))) > abs(corr(all_time(idx)', V2(:, i2(2))))
                t_dom(sim_num, i1, k) = 1;
            end
                    figure(f1);
                    subplot(nlambda, ntime, fig_idx)
                    scatter(V2(:, i2(1)), V2(:, i2(2)),50, all_p(idx), '.')
                    title(sprintf('t=%d, \\lambda=%4.0f', i1, lambda))
                    %         axis off
            
                    figure(f2);
                    subplot(nlambda, ntime, fig_idx)
                    scatter(V2(:, i2(1)), V2(:, i2(2)),50, all_time(idx), '.')
                    title(sprintf('t=%d, \\lambda=%4.0f', i1, lambda))
                    %         axis off
            
                    fig_idx = fig_idx + 1;
            
        end
        
    end
end


%%

C = mean(sqrt(log(eval1)./log(eval2)), 3);
C2 = repmat(t_range_all, nlambda, 1).*N./repmat(lambda_all', 1, ntime);

% X = unique(lam_vec);
% Y = unique(t_vec);
% C = zeros(length(X), length(Y));
% C2 = zeros(length(X), length(Y));
%
% for i=1:length(eval1_vec)
%     i1 = find(X == lam_vec(i));
%     i2 = find(Y == t_vec(i));
%     C(i1, i2) = C(i1, i2) + sqrt(log(eval1_vec(i))./log(eval2_vec(i)));
%     C2(i1, i2) = C2(i1, i2) + max(lam_vec(i)*t_vec(i), 1/(lam_vec(i)*t_vec(i)));
% end
%
% C  = C / niter;
% C2  = C2 / niter;


figure;
imagesc(lambda_all, t_range_all, C')
set(gca, 'xdir','normal')
set(gca, 'ydir','normal')
colorbar


figure;
imagesc(lambda_all, t_range_all, (C2'))
set(gca, 'xdir','normal')
set(gca, 'ydir','normal')
colorbar

figure;
imagesc(lambda_all, t_range_all, t_dom')
set(gca, 'xdir','normal')
set(gca, 'ydir','normal')
colorbar