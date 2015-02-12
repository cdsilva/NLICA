clear all
close all

rng(123);
markersize = 500;


%% simulations

lambda_all = [100 2500 6400];
s_all = sqrt(lambda_all);
comp2 = [3 4 5];

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
%nbins = 128;
nbins = 32;
hist_all = zeros(nbins, nsteps*nsim);
hist_left = zeros(nbins, nsteps*nsim);
hist_right = zeros(nbins, nsteps*nsim);
all_time = zeros(1, nsteps*nsim);
all_p = zeros(1, nsteps*nsim);


for sim_num = 1:length(lambda_all)
    
    lambda = lambda_all(sim_num);
    s = s_all(sim_num);
    
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
    
    
    if corr(V2(:,2), all_time(idx)') < 0
        V2(:,2) = -V2(:,2);
    end
    if corr(V2(:,3), all_p(idx)') < 0
        V2(:,3) = -V2(:,3);
    end
    
    %%
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V2, eps_med_scale);
    
    make_fig(3,2);
    make_colored_bars(diag(D2(2:end,2:end)), res(2:end))
    if sim_num ~= length(lambda_all)
        colorbar off
    end
    xlabel('k')
    ylabel('\mu_k')
    axis square
    print(sprintf('chemotaxis%d_evals.eps', sim_num),'-depsc')
    
    make_fig(3,2);
    scatter(V2(:,2),V2(:,3),50, all_p(idx), '.')
    if sim_num == length(lambda_all)
        cbar = colorbar('peer',gca);
        set(get(cbar,'xlabel'),'String','p');
    end
    xlabel('\phi_1')
    ylabel('\phi_2')
    axis tight
    axis square
    print(sprintf('chemotaxis%d_embed_bad.eps', sim_num),'-depsc')
    
    make_fig(3,2);
    scatter(V2(:,2),V2(:,comp2(sim_num)),50, all_p(idx), '.')
    if sim_num == length(lambda_all)
        cbar = colorbar('peer',gca);
        set(get(cbar,'xlabel'),'String','p');
    end
    xlabel('\phi_1')
    ylabel(sprintf('\\phi_%d', comp2(sim_num)-1))
    axis tight
    axis square
    print(sprintf('chemotaxis%d_embed_good.eps', sim_num),'-depsc')
    
    %
    % figure;
    % bar(res(2:end))
    % ylabel('cross-validation error')
    % % print(sprintf('lle_errors%d', sim_num),'-depsc')
    %
    % figure;
    % bar(diag(D2(2:end,2:end)))
    % ylabel('\lambda_k')
    % % print(sprintf('lle_evals%d', sim_num),'-depsc')
    %
    % figure;
    % scatter(V2(:,2),V2(:,3),markersize,'b', '.')
    % xlabel('$\phi_1$', 'interpreter','latex')
    % ylabel('$\phi_2$', 'interpreter','latex')
    % % print(sprintf('lle_embed1_%d', sim_num),'-depsc')
    %
    % figure;
    % scatter(V2(:,2),V2(:,4),markersize,'b', '.')
    % xlabel('$\phi_1$', 'interpreter','latex')
    % ylabel('$\phi_3$', 'interpreter','latex')
    % % print(sprintf('lle_embed2_%d', sim_num),'-depsc')
    
end


%%


% lambda_all = [100 2500 6400];
lambda_all = linspace(1, 4000, 40);
s_all = sqrt(lambda_all);

eval1 = zeros(size(lambda_all));
eval2 = zeros(size(lambda_all));

for sim_num = 1:length(lambda_all)
    
    lambda = lambda_all(sim_num);
    s = s_all(sim_num);
    
    x_hist = linspace(min_start-s*tmax, max_start+s*tmax, nbins);
    
    
    %
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
    
    
    % use data after initial "relaxation" time
    idx = (all_time > 1);
    
    % compute EMD between histograms
    
    W2 = zeros(nsteps*nsim);
    
    D_for_emd = pdist2(x_hist', x_hist');
    
    for i1=1:nsteps*nsim
        for i2=1:i1-1
            
            emd_dist = sum(abs(cumsum(hist_all(:,i1)) - cumsum(hist_all(:,i2))));
            
            W2(i1, i2) = emd_dist^2;
            W2(i2, i1) = W2(i1, i2);
        end
    end
    
    % dmaps on EMD
    W2 = W2(idx,idx);
    eps2 = median(W2(:));
    [V2, D2] = dmaps(W2, eps2, 5);
    
    
    % compute regression
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V2, eps_med_scale);
    
    [~, idx] = sort(res(2:end), 'descend');
    idx = idx + 1;
    if idx(1) > idx(2)
        tmp = idx(1);
        idx(1) = idx(2);
        idx(2) = tmp;
    end
    
    eval1(sim_num) = D2(idx(1), idx(1));
    eval2(sim_num) = D2(idx(2), idx(2));
    
    
    
end

%%
% figure;
% plotyy(lambda_all, log(eval1)./log(eval2), lambda_all, log(lambda_all.*s_all* tmax))

L = s_all.*tmax;
tau_run = 1./lambda_all;
tau_drift = L./s_all;
tau_diff = L.^2.*lambda_all./(s_all.^2);

make_fig(4,3);
[ax,h1,h2] = plotyy(lambda_all, log(eval1)./log(eval2), lambda_all, tau_drift./tau_diff, 'plot','semilogy');
set(h1, 'marker', '.', 'color','b')
set(h2, 'linestyle', '-', 'color','r')
% set(h2, 'yscale', 'log')
set(ax(1), 'ycolor','b')
set(ax(2), 'ycolor','r')
ylabel(ax(1), 'log(\mu_{i_1})/log(\mu_{i_2})');
ylab = ylabel(ax(2), '\tau_{diff}/\tau_{drift}');
set(ylab, 'Units', 'Normalized', 'Position', [1.02, 0.5, 0]);
xlabel('\lambda')
% axis tight
% axis square
print('chemotaxis_compare_timescales_evals.eps','-depsc')




