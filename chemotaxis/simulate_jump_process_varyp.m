clear all
close all

%% simulations

% rate of switching velocity
nsiml = 10;
lambda = logspace(1,3,nsiml);
%lambda = linspace(10,1000,nsiml);
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
nsimp = 10;
pinit = linspace(0.1, 0.9, nsimp); 

% histogram parameters
nbins = 128;
x_hist = linspace(-s*tmax, s*tmax, nbins);
hist_all = zeros(nbins, nsteps*nsimp*nsiml);
hist_left = zeros(nbins, nsteps*nsimp*nsiml);
hist_right = zeros(nbins, nsteps*nsimp*nsiml);
all_time = zeros(1, nsteps*nsimp*nsiml);
all_p = zeros(1, nsteps*nsimp*nsiml);
all_lambda = zeros(1, nsteps*nsimp*nsiml);

subplot_dim1 = floor(sqrt(nsimp));
subplot_dim2 = ceil(nsimp / subplot_dim1);

figure;
for k=1:nsiml
    for j=1:nsimp
        % simulate velocity jump process
        [pos, vel, time] = simulate_jumps(N, tmax, dt, lambda(k), s, pinit(j));

        % histogram of all particle positions as a function of time
        for i=1:nsteps

            [hist_all(:,(k-1)*nsteps*nsimp+(j-1)*nsteps+i), ~] = hist(pos(:,i),x_hist);

            idx = (vel(:,i) == -1);
            [hist_left(:,(k-1)*nsteps*nsimp+(j-1)*nsteps+i), ~] = hist(pos(idx,i),x_hist);

            idx = (vel(:,i) == 1);
            [hist_right(:,(k-1)*nsteps*nsimp+(j-1)*nsteps+i), ~] = hist(pos(idx,i),x_hist);

        end

        % store times and probabilities
        all_time((k-1)*nsteps*nsimp+(j-1)*nsteps+1:(k-1)*nsteps*nsimp+j*nsteps) = time;
        all_p((k-1)*nsteps*nsimp+(j-1)*nsteps+1:(k-1)*nsteps*nsimp+j*nsteps) = pinit(j);
        all_lambda((k-1)*nsteps*nsimp+(j-1)*nsteps+1:(k-1)*nsteps*nsimp+j*nsteps) = lambda(k);
        
        % plot histograms
        %subplot(subplot_dim1, subplot_dim2, j)
%         figure;
%         plot(x_hist,hist_all(:,(j-1)*nsteps+1:j*nsteps))
%         title(sprintf('p = %0.2f',pinit(j)))

    end
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

figure;
scatter(V(:,2),V(:,3),50,log(all_lambda(idx)), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('histograms: colored by \lambda')

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

for i=1:nsteps*nsimp*nsiml
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
eps2 = 10*median(W2(:));

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

figure;
scatter(V2(:,2),V2(:,3),50,log(all_lambda(idx)), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('scattering transform: colored by \lambda')

figure;
scatter3(V2(:,2),V2(:,3),V2(:,5),50,all_time(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('scattering transform: colored by time')

figure;
scatter3(V2(:,2),V2(:,3),V2(:,5),50,all_p(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('scattering transform: colored by p')

figure;
scatter3(V2(:,2),V2(:,3),V2(:,5),50,log(all_lambda(idx)), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('scattering transform: colored by \lambda')

figure;
scatter3(all_time(idx),all_p(idx),all_lambda(idx),500,V2(:,2),'.')
xlabel('time')
ylabel('p')
zlabel('\lambda')
title('colored by \phi_2')

figure;
scatter3(all_time(idx),all_p(idx),all_lambda(idx),500,V2(:,3),'.')
xlabel('time')
ylabel('p')
zlabel('\lambda')
title('colored by \phi_3')

%% compute wavelet transform for histograms

% % order of wavelet transform
% wavelet_N = 4;
% 
% W2 = zeros(nsteps*nsimp);
% 
% for i1=1:nsteps*nsimp
%     for i2=1:i1-1
%         [coeff,L] = wavedec(hist_all(:,i1)-hist_all(:,i2),wavelet_N,'db12');
%         L2 = zeros(L(1),1);
%         for j=2:length(L)-1
%             L2 = [L2; 2^-((j-2)*(1+1/2))*ones(L(j),1)];
%         end
%         
%         W2(i1, i2) = sum(L2.*abs(coeff));
%         W2(i2, i1) = W2(i1, i2);
%     end
% end
% 
% %% dmaps on wavelet coefficients
% 
% W2 = W2(idx,idx);
% 
% eps2 = median(W2(:));
% 
% [V2, D2] = dmaps(W2, eps2, 10);
% 
% figure;
% scatter(V2(:,2),V2(:,3),50,all_time(idx), '.')
% xlabel('\phi_2')
% ylabel('\phi_3')
% title('wavelet transform: colored by time')
% 
% figure;
% scatter(V2(:,2),V2(:,3),50,all_p(idx), '.')
% xlabel('\phi_2')
% ylabel('\phi_3')
% title('wavelet transform: colored by p')

%% compute EMD for histograms

addpath('FastEMD');

W2 = zeros(nsteps*nsimp*nsiml);
D_for_emd = pdist2(x_hist', x_hist');

for i1=1:nsteps*nsimp*nsiml
    for i2=1:i1-1
        [emd_dist, ~]= emd_hat_gd_metric_mex(hist_all(:,i1), hist_all(:,i2), D_for_emd);
        
        
        W2(i1, i2) = emd_dist^2;
        W2(i2, i1) = W2(i1, i2);
    end
end

%% dmaps on EMD

eps2 = median(W2(:));

[V2, D2] = dmaps(W2, eps2, 10);

figure;
scatter(V2(:,2),V2(:,3),50,all_time(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('EMD transform: colored by time')

figure;
scatter(V2(:,2),V2(:,3),50,all_p(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('EMD transform: colored by p')

figure;
scatter(V2(:,2),V2(:,3),50,log(all_lambda(idx)), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('EMD transform: colored by \lambda')

figure;
scatter3(V2(:,2),V2(:,3),V2(:,5),50,all_time(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('EMD transform: colored by time')

figure;
scatter3(V2(:,2),V2(:,3),V2(:,5),50,all_p(idx), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('EMD transform: colored by p')

figure;
scatter3(V2(:,2),V2(:,3),V2(:,5),50,log(all_lambda(idx)), '.')
xlabel('\phi_2')
ylabel('\phi_3')
title('EMD transform: colored by \lambda')

figure;
scatter3(all_time(idx),all_p(idx),all_lambda(idx),500,V2(:,2),'.')
xlabel('time')
ylabel('p')
zlabel('\lambda')
title('colored by \phi_2')

figure;
scatter3(all_time(idx),all_p(idx),all_lambda(idx),500,V2(:,3),'.')
xlabel('time')
ylabel('p')
zlabel('\lambda')
title('colored by \phi_3')


