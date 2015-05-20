% example of using local linear regression to find unique eigendirections
% for a two-dimensional strip data set

%% define initial parameters

% length along x-direction
Lx = 4;

% length along y-direction
Ly = 1;

% number of data points
n = 2000;

%% generate data

% intialize random number generator
rng(321);

% data uniformly distributed on Lx x Ly strip
data = rand(n, 2);
data(:,1) = data(:,1) * Lx;
data(:,2) = data(:,2) * Ly;

%% diffusion maps

% pairwise distance matrix
W = squareform(pdist(data));

% kernel scale
eps = 0.15*sqrt(2);

% compute embedding coordinates
[V, D] = dmaps(W, eps, 10);

%% local linear regression

% regression kernel scale
eps_med_scale = 3;

% compute cross-validation error in fitting diffusion maps eigenvectors 
% as a function of previous eigenvectors
res = compute_residuals_DMAPS(V, eps_med_scale);

%% make plots

% plot original data
figure;
scatter(data(:,1),data(:,2),25, 'b','.')
axis equal
xlabel('z_1')
ylabel('z_2')

% plot diffusion maps eigenvalues, colored by cross-validation error
figure;
colored_bars(diag(D(2:end,2:end)), res(2:end))
xlabel('k')
ylabel('\mu_k')

