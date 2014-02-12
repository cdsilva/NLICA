% this code simulates the data from "Non LInear INdependent Component
% Analysis"

clear all
close all

%% define parameters
N = 2000;
dim = 2;
Nc = 200;
dt = 0.01;

%% do simulations

x = rand(N, dim);

f = @(x) [x(:,1)+x(:,2).^3 x(:,2)-x(:,1).^3];
y = f(x);

%% plot results

% plot variables, color by time
figure;
plot(x(:,1),x(:,2),'.')

figure;
plot(y(:,1),y(:,2),'.')

%% calculate clouds

xclouds = dt * randn(Nc, dim, N);
for i=1:N
    xclouds(:,:,i) = repmat(x(i,:), Nc, 1) + xclouds(:,:,i);
end
ind = find(xclouds > 1);
xclouds(ind) = 2 - xclouds(ind);
ind = find(xclouds < 0);
xclouds(ind) = -xclouds(ind);

yclouds = zeros(size(xclouds));
for i=1:N
    yclouds(:,:,i) = f(xclouds(:,:,i));
end

%% calculate inverse covariances
[inv_c, ~, ranks] = covariances2(yclouds, dim);

%% NIV
if size(y,1) > 4000
    disp('Too much data; NIV will fail')
    return
end
neigs = 10;
[V, D] = NIV(y, inv_c, 0, neigs, 0);

figure;
scatter(x(:,1),x(:,2),200,V(:,2),'.')

figure;
scatter(x(:,1),x(:,2),200,V(:,3),'.')

%% dmaps
W = squareform(pdist(y)).^2;
[V, D] = dmaps(W, median(W(:)), neigs);

figure;
scatter(x(:,1),x(:,2),200,V(:,2),'.')

figure;
scatter(x(:,1),x(:,2),200,V(:,3),'.')

%% do simulation with noise

dim = 3;

x = zeros(N, dim);
x(:, 1:2) = rand(N, 2);

f = @(x) [x(:,1)+x(:,2).^3 x(:,2)-x(:,1).^3 x(:,3)];

y = f(x);

%% do simulations

xclouds = zeros(Nc, dim, N);

drift = @(t, x) [0; 0; (1-x(3))];
diffn = @(t, x) [1 0 0; 0 1 0; 0 0 100];
    
dt = dt/Nc;

for i=1:N
    SDE = sde(drift, diffn, 'StartState', x(i,:));
    [xclouds(:,:,i), ~, ~] = SDE.simulate(Nc, 'DeltaTime', dt);
end


