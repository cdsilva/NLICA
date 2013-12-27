clear all
close all

%% simulate 2 wells

% simulation parameters
N = 4000;
% time step
dt = 0.02;
% step intervals at which to record data
step_size = 5;
% dimension
dim = 1;
% diffusivity
D = 0.4;
% initial condition
x0 = 0.3;
% integrate SDE
[data, x, g] = simulate_SDE2(2*step_size*N, dt, D, x0);
% subsample data
data = data(1:step_size:end);

% construct two data sets from simulation
noverlap = 10;
data1 = data(1:N);
data2 = data(N-noverlap:end);
N1 = length(data1);
N2 = length(data2);

% create function of data
f = @(x) x.^3-2*x;  
data2b = f(data2);

%% plot
figure;
plot(x, g, '.')
hold on
plot(data1, (1:length(data1))/length(data1),'.','color','k')
title('Sampled data 1')

figure;
plot(x, g, '.')
hold on
plot(data2, (1:length(data2))/length(data2),'.','color','k')
title('Sampled data 2')

figure;
[N, X] = hist(data1, 50);
plot(X, N / sum(N))
hold on
[N, X] = hist(data2, 50);
plot(X, N / sum(N), 'r')
legend('data 1','data 2','location','best')
title('Histogram of data')

%% dmaps

evec_idx = 2;

W = squareform(pdist(data1)).^2;
%eps = 0.5;
eps = median(median(W));
alpha = 1;
[V1, D1] = dmaps_weight(W, alpha, eps, 5);
V1(:, evec_idx) = V1(:,evec_idx) * sign(corr(V1(:,evec_idx), data1));

W = squareform(pdist(data2)).^2;
%eps = 0.1;
eps = median(median(W));
[V2, D2] = dmaps_weight(W, alpha, eps, 5);
V2(:, evec_idx) = V2(:,evec_idx) * sign(corr(V2(:,evec_idx), data2));

W = squareform(pdist(data2b)).^2;
%eps = 0.1;
eps = median(median(W));
[V2b, D2b] = dmaps_weight(W, alpha, eps, 5);
V2b(:, evec_idx) = V2b(:,evec_idx) * sign(corr(V2b(:,evec_idx), data2));

figure;
plot(data1, V1(:,evec_idx)/length(data1), '.')
hold on
plot(data2, V2(:,evec_idx)/length(data2), '.r')
plot(data2, V2b(:,evec_idx)/length(data2), '.g')
xlabel('x')
ylabel('\phi_2')
legend('data 1','data 2','data 2 (function)','location','best')
title('DMAPS embedding')

%% align embeddings

% P = polyfit(V1(end-noverlap:end,evec_idx),V2(1:noverlap+1,evec_idx),1);
% 
% figure; 
% plot(V1(end-noverlap:end,evec_idx),V2(1:noverlap+1,evec_idx),'.')
% xlabel('\phi_2 (data 1)')
% ylabel('\phi_2 (data 2)')
% 
% figure;
% plot(data1, polyval(P, V1(:,2)),'.')
% hold on
% plot(data2,V2(:,2),'.r')
% xlabel('x')
% ylabel('\phi_2')

%% calculate covariances

inv_c1 = zeros(dim, dim, N1);
inv_c2 = zeros(dim, dim, N2);
inv_c2b = zeros(dim, dim, N2);

ncov_points = 10;
for i=1:N1
    Z = data1(max(1, i-ncov_points):min(N1,i+ncov_points));
    inv_c1(:,:,i) = 1/var(Z);
end
for i=1:N2
    Z = data2(max(1, i-ncov_points):min(N2,i+ncov_points));
    inv_c2(:,:,i) = 1/var(Z);
    Z = data2b(max(1, i-ncov_points):min(N2,i+ncov_points));
    inv_c2b(:,:,i) = 1/var(Z);
end

%% NLICA
alpha = 1;
neigs = 5;

[V1, D1] =  NLICA2(data1, inv_c1, neigs, alpha); 
V1(:, evec_idx) = V1(:,evec_idx) * sign(corr(V1(:,evec_idx), data1));

[V2, D2] =  NLICA2(data2, inv_c2, neigs, alpha); 
V2(:, evec_idx) = V2(:,evec_idx) * sign(corr(V2(:,evec_idx), data2));

[V2b, D2b] =  NLICA2(data2b, inv_c2b, neigs, alpha); 
V2b(:, evec_idx) = V2b(:,evec_idx) * sign(corr(V2b(:,evec_idx), data2));

figure;
plot(data1, V1(:,evec_idx)/length(data1), '.')
hold on
plot(data2, V2(:,evec_idx)/length(data2), '.r')
plot(data2, V2b(:,evec_idx)/length(data2), '.g')
xlabel('x')
ylabel('\phi_2')
legend('data 1','data 2','data 2 (function)','location','best')
title('NIV embedding')

%% variance of NIV

var_niv1 = zeros(N1,1);
for i=1:N1
    var_niv1(i) = var(V1(max(1, i-ncov_points):min(N1,i+ncov_points),evec_idx));
end

figure;
plot(V1(:, evec_idx),var_niv1, '.')

figure; 
plot(var_niv1, '.')

[~, I] = sort(V1(:, evec_idx));
figure;
plot(V1(I, evec_idx),smooth(var_niv1(I), 25), '.')

%% histogram

figure;
hist(V1(:,evec_idx),300)

figure;
plot(cumsum(