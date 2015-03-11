clear all
close all

rng(123);


%% simulations

% nlambda = 10;
% lambda_all = linspace(1, 400, nlambda);
lambda_all = 1:10:200;
nlambda = length(lambda_all);

% number of particles
N = 1000;

% ntmax = 10;
% tmax_all = linspace(10, 400, ntmax);
tmax_all = 10:10:200;
ntmax = length(tmax_all);

%initial probility of a particle moving to the right (vel=+1)
nsim = 10;
pinit = linspace(0.1, 0.9, nsim);

% histogram parameters
nbins = 32;

niter = 3;

eval1 = zeros(nlambda, ntmax, niter);
eval2 = zeros(nlambda, ntmax, niter);
t_dom = zeros(nlambda, ntmax, niter);

% f1 = figure;
% f2 = figure;
% fig_idx = 1;


for j1=1:nlambda
    j1
    for j2=1:ntmax
        for k=1:niter
            
            % rate of switching velocity
            lambda = lambda_all(j1);
            % speed of particles
            s = sqrt(lambda);
            
            % maximum simulation time
            tmax = tmax_all(j2);
            %             tmin = 0.1*tmax;
            tmin = 0;
            
            % time step
            dt = tmax / 10;
            nsteps = floor(tmax / dt+1e-6) + 1;
            
            
            % histogram parameters
            x_hist = linspace(-s*tmax, s*tmax, nbins);
            hist_all = zeros(nbins, nsteps*nsim);
            hist_left = zeros(nbins, nsteps*nsim);
            hist_right = zeros(nbins, nsteps*nsim);
            all_time = zeros(1, nsteps*nsim);
            all_p = zeros(1, nsteps*nsim);
            
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
            idx = (all_time > tmin);
            
            % compute EMD between histograms
            
            W2 = zeros(nsteps*nsim);
            
            D_for_emd = pdist2(x_hist', x_hist');
            
            for i1=1:nsteps*nsim
                for idx2=1:i1-1
                    
                    emd_dist = sum(abs(cumsum(hist_all(:,i1)) - cumsum(hist_all(:,idx2))));
                    
                    W2(i1, idx2) = emd_dist^2;
                    W2(idx2, i1) = W2(i1, idx2);
                end
            end
            
            
            % dmaps on EMD
            
            W2 = W2(idx,idx);
            eps2 = median(W2(:));
            [V2, D2] = dmaps(W2, eps2, 5);
            
            
            %
            eps_med_scale = 3;
            res = compute_residuals_DMAPS(V2, eps_med_scale);
            
            %         i2 = find(res > 0.3);
            [~, idx2] = sort(res, 'descend');
            idx2 = sort(idx2(1:2));
            %         i2 = [2 3];
            
            eval1(j1, j2, k) = D2(idx2(1), idx2(1));
            eval2(j1, j2, k) = D2(idx2(2), idx2(2));
            
            if abs(corr(all_time(idx)', V2(:, idx2(1)))) > abs(corr(all_time(idx)', V2(:, idx2(2))))
                t_dom(j1, j2, k) = 1;
            end
            %
            %             figure(f1);
            %             subplot(nlambda, ntmax, fig_idx)
            %             scatter(V2(:, idx2(1)), V2(:, idx2(2)),50, all_p(idx), '.')
            %             title(sprintf('t=%d, \\lambda=%4.0f', tmax, lambda))
            %             %         axis off
            %
            %             figure(f2);
            %             subplot(nlambda, ntmax, fig_idx)
            %             scatter(V2(:, idx2(1)), V2(:, idx2(2)),50, all_time(idx), '.')
            %             title(sprintf('t=%d, \\lambda=%4.0f', tmax, lambda))
            %             %         axis off
            %
            %             fig_idx = fig_idx + 1;
            
        end
    end
end


%%

if niter > 1
    C = mean(sqrt(log(eval1)./log(eval2)), 3);
else
    C = sqrt(log(eval1)./log(eval2));
end

make_fig(4,3);
imagesc(lambda_all, tmax_all, C')
set(gca, 'xdir','normal')
set(gca, 'ydir','normal')
h = colorbar;
set(get(h,'ylabel'),'String', '$\sqrt{\log \mu_{i_1} / \log \mu_{i_2}}$', 'interpreter','latex');
hold on
plot(1:max(lambda_all), N./(1:max(lambda_all)), '-w', 'linewidth',2)
xlabel('\lambda')
ylabel('t_{max}')
saveas(gcf, 'tmax_lambda_transition', 'epsc')

figure;
if niter > 1
    imagesc(lambda_all, tmax_all, mean(t_dom, 3)')
else
    imagesc(lambda_all, tmax_all, t_dom')
end
set(gca, 'xdir','normal')
set(gca, 'ydir','normal')
colorbar
