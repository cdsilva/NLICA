clear all
close all

rng(123);
markersize = 500;


%% simulations

% rate of switching velocity
lambda = 100;
% speed of particles
s = 10;
sim_num = 1;

% % rate of switching velocity
% lambda = 2500;
% % speed of particles
% s = 50;
% sim_num = 2;

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

W = squareform(pdist(hist_all(:,idx)')).^2;
eps = median(W(:));

[V, D] = dmaps(W, eps, 10);

figure;
scatter(V(:,2),V(:,3),markersize,all_time(idx), '.')
xlabel('$\phi_1$', 'interpreter','latex', 'fontsize', 20)
ylabel('$\phi_2$', 'interpreter','latex', 'fontsize', 20)
%title('histograms: colored by time')
h = colorbar;
set(get(h,'xlabel'),'String', 't', 'fontsize', 20);
%print(sprintf('rawhist_t_%d', lambda), '-r300','-djpeg')
% saveas(gcf, sprintf('rawhist_t_%d', lambda), 'epsc')

figure;
scatter(V(:,2),V(:,3),markersize,all_p(idx), '.')
xlabel('$\phi_1$', 'interpreter','latex', 'fontsize', 20)
ylabel('$\phi_2$', 'interpreter','latex', 'fontsize', 20)
%title('histograms: colored by p')
h = colorbar;
set(get(h,'xlabel'),'String', 'p', 'fontsize', 20);
%print(sprintf('rawhist_p_%d', lambda), '-r300','-djpeg')
% saveas(gcf, sprintf('rawhist_p_%d', lambda), 'epsc')

corr(V(:,2),all_time(idx)')

corr(V(:,3),all_p(idx)')

%% compute EMD between histograms

addpath('FastEMD');

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

rmpath('FastEMD');

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

%%
eps_med_scale = 3;
res = compute_residuals_DMAPS(V2, eps_med_scale);

figure;
bar(res(2:end))
ylabel('cross-validation error')
print(sprintf('lle_errors%d', sim_num),'-depsc')

figure;
bar(diag(D2(2:end,2:end)))
ylabel('\lambda_k')
print(sprintf('lle_evals%d', sim_num),'-depsc')

figure;
scatter(V2(:,2),V2(:,3),markersize,'b', '.')
xlabel('$\phi_1$', 'interpreter','latex')
ylabel('$\phi_2$', 'interpreter','latex')
print(sprintf('lle_embed1_%d', sim_num),'-depsc')

figure;
scatter(V2(:,2),V2(:,4),markersize,'b', '.')
xlabel('$\phi_1$', 'interpreter','latex')
ylabel('$\phi_3$', 'interpreter','latex')
print(sprintf('lle_embed2_%d', sim_num),'-depsc')

return
%%

num_hist = 6;

figure;
set(gcf,'PaperPositionMode','auto');
scatter(V2(:,2),V2(:,3),markersize,all_time(idx), '.');
ax = gca;
axis_lim = [-0.4 0.25 -0.35 0.4];
axis(axis_lim)
xlabel('$\phi_1$', 'interpreter','latex', 'fontsize', 20)
ylabel('$\phi_2$', 'interpreter','latex', 'fontsize', 20)
h = colorbar;
set(get(h,'xlabel'),'String', 't', 'fontsize', 20);
curr_ax = get(ax, 'position');
curr_ax(1) = curr_ax(1) + 0.02;
set(ax, 'position', curr_ax);
hist_all2 = hist_all(:, idx);
hold on
i = zeros(num_hist, 1);
i(1) = find(all_time(idx) == min(all_time(idx)) & all_p(idx) == min(all_p(idx)));
i(2) = find(all_time(idx) == max(all_time(idx)) & all_p(idx) == min(all_p(idx)));
i(3) = find(all_time(idx) == min(all_time(idx)) & all_p(idx) == max(all_p(idx)));
i(4) = find(all_time(idx) == max(all_time(idx)) & all_p(idx) == max(all_p(idx)));
i(5) = find(all_time(idx) == floor(mean(all_time(idx))-2) & all_p(idx) == min(all_p(idx)));
i(6) = find(all_time(idx) == floor(mean(all_time(idx))-2) & all_p(idx) == max(all_p(idx)));

points_x = curr_ax(1) + (V2(i, 2) - axis_lim(1)) / (axis_lim(2)-axis_lim(1)) * curr_ax(3);
points_y = curr_ax(2) + (V2(i, 3) - axis_lim(3)) / (axis_lim(4)-axis_lim(3)) * curr_ax(4);

ax_size = 0.15;
ax_pos = [.17 .16 ax_size ax_size;
    .6 .16 ax_size ax_size;
    .17 .75 ax_size ax_size;
    .6 .75 ax_size ax_size;
    .375 .16 ax_size ax_size;
    .375 .75 ax_size ax_size];

arrow_pos = ax_pos(:, 1:2);
arrow_pos(1:4,1) = arrow_pos(1:4,1) + ax_size/2;
arrow_pos(1:2,2) = arrow_pos(1:2,2) + ax_size;
arrow_pos(5:6,1) = arrow_pos(5:6,1) + ax_size/2;
arrow_pos(5,2) = arrow_pos(5,2) + ax_size;

xx = linspace(min_start-s*tmax, max_start+s*tmax, 500);

for j=1:num_hist
    axes('position',ax_pos(j, :));
    yy = interp1(x_hist, hist_all2(:, i(j)),xx, 'pchip');
    plot(xx, yy)
    %bar(x_hist, hist_all2(:, i(j)))
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
end
for j=1:num_hist
    annotation('arrow', [points_x(j) arrow_pos(j, 1)], [points_y(j) arrow_pos(j, 2)])
end
%print(sprintf('EMD_withhist_t_%d', lambda), '-r300','-djpeg')
saveas(gcf, sprintf('EMD_withhist_t_%d', lambda), 'epsc')


figure;
set(gcf,'PaperPositionMode','auto');
scatter(V2(:,2),V2(:,3),markersize,all_p(idx), '.');
ax = gca;
axis(axis_lim)
xlabel('$\phi_1$', 'interpreter','latex', 'fontsize', 20)
ylabel('$\phi_2$', 'interpreter','latex', 'fontsize', 20)
h = colorbar;
set(get(h,'xlabel'),'String', 'p', 'fontsize', 20);
curr_ax = get(ax, 'position');
curr_ax(1) = curr_ax(1) + 0.02;
set(ax, 'position', curr_ax);
hold on

for j=1:num_hist
    axes('position',ax_pos(j, :));
    yy = interp1(x_hist, hist_all2(:, i(j)),xx, 'pchip');
    plot(xx, yy)
    %bar(x_hist, hist_all2(:, i(j)))
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
end
for j=1:num_hist
    annotation('arrow', [points_x(j) arrow_pos(j, 1)], [points_y(j) arrow_pos(j, 2)])
end
%print(sprintf('EMD_withhist_p_%d', lambda), '-r300','-djpeg')
saveas(gcf, sprintf('EMD_withhist_p_%d', lambda), 'epsc')

