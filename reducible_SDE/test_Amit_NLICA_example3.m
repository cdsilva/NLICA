% this code simulates the data from "Non LInear INdependent Component
% Analysis"

clear all
close all

res = '-r300';
fmt = '-djpeg';

rng(321);

%% define parameters

% number of data points
N = 3500;

% dimension
dim = 2;

% number of points in each cloud
Nc = 5000;

% dt for clouds
dt1 = 0.01;

%% do simulations

x = rand(N, dim);

alpha = 20;
sin_scale = 0.2;
f = @(x) [x(:,1)+x(:,2).^3 x(:,2)-x(:,1).^3 sin_scale*sin(alpha*(x(:,1)+x(:,2).^3))+sin_scale*sin(alpha*(x(:,2)-x(:,1).^3))];
%f = @(x) [x(:,1) x(:,2) 0.1*sin(alpha*x(:,1))+0.1*sin(alpha*x(:,2))];
y = f(x);

%% plot results

figure;
plot(x(:,1),x(:,2),'.')
xlabel('x_1')
ylabel('x_2')
%print('xdata', fmt, res)

figure;
plot3(y(:,1),y(:,2), y(:, 3),'.')
xlabel('y_1')
ylabel('y_2')
%print('ydata', fmt, res)

%% calculate clouds

xclouds = sqrt(dt1) * randn(Nc, dim, N);
for i=1:N
    xclouds(:,:,i) = repmat(x(i,:), Nc, 1) + xclouds(:,:,i);
end
ind = find(xclouds > 1);
xclouds(ind) = 2 - xclouds(ind);
ind = find(xclouds < 0);
xclouds(ind) = -xclouds(ind);

yclouds = zeros(Nc, 3, N);
for i=1:N
    yclouds(:,:,i) = f(xclouds(:,:,i));
end

%% calculate inverse covariances
[inv_c, ~, ranks] = covariances2(yclouds, dim);
inv_c = inv_c * dt1;

%% NIV
if size(y,1) > 4000
    disp('Too much data; NIV will fail')
    return
end
neigs = 10;
eps = 2*0.005;
%eps = 0;
[V, D, ~] = NIV(y, inv_c, eps, neigs, 0);
V = decouple_dmaps(V);

figure;
scatter(x(:,1),x(:,2),200,V(:,2),'.')
xlabel('x_1')
ylabel('x_2')
%print('xdata_colored_NIV1', fmt, res)

figure;
scatter(x(:,1),x(:,2),200,V(:,3),'.')
xlabel('x_1')
ylabel('x_2')
%print('xdata_colored_NIV2', fmt, res)

figure;
plot(V(:,2),V(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')
%print('embedding_niv', fmt, res)

figure;
bar(-2*log(diag(D))/(pi^2 * eps/2))

 
%% dmaps
W = squareform(pdist(y)).^2;
[V, D] = dmaps(W, eps, neigs);
V = decouple_dmaps(V);

figure;
scatter(x(:,1),x(:,2),200,V(:,2),'.')
%print('xdata_colored_DMAPS1', fmt, res)

figure;
scatter(x(:,1),x(:,2),200,V(:,3),'.')
%print('xdata_colored_DMAPS2', fmt, res)

figure;
plot(V(:,2),V(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')
%print('embedding_dmaps', fmt, res)

figure;
bar(-2*log(diag(D))/(pi^2 * eps/2))


%% do simulation with noise

dt2 = dt1 * 2;

z = y;
dim = 3;

%zclouds = dt2 * randn(Nc, dim, N);
zclouds = dt2 * (2*rand(Nc, dim, N)-1);
for i=1:N
    zclouds(:,:,i) = zclouds(:,:,i) + yclouds(:,:,i);
end

% calculate inverse covariances
[inv_c, ~, ranks] = covariances2(zclouds, dim);
inv_c = inv_c * dt1;

% NIV
if size(z,1) > 4000
    disp('Too much data; NIV will fail')
    return
end

[V, D, ~] = NIV(z, inv_c, eps, neigs, 0);
V = decouple_dmaps(V);

figure;
scatter(x(:,1),x(:,2),200,V(:,2),'.')
xlabel('x_1')
ylabel('x_2')
%print('xdata_noise2_colored_NIV1', fmt, res)

figure;
scatter(x(:,1),x(:,2),200,V(:,3),'.')
xlabel('x_1')
ylabel('x_2')
%print('xdata_noise2_colored_NIV2', fmt, res)

figure;
plot(V(:,2),V(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')
%print('embedding_noise2', fmt, res)

figure;
bar(-2*log(diag(D))/(pi^2 * eps/2))

%%

zclouds2 = zeros(Nc*N, dim);
for i=1:N
    zclouds2((i-1)*Nc+1:i*Nc,:) = zclouds(:, :, i);
end

%%
ind = 1885;

nscales = 5;
eps_nn = sqrt(3)*linspace(dt2*2, dt2/4, nscales);

%nvecs = 7;
%proj_vec = randn(3, nvecs);
nvecs = 3;
proj_vec = eye(3);

colors = 'brgky';

D_all = zeros(nscales, nvecs);

%%
figure;
for i=1:nscales
    idx = find(sum((zclouds2 - repmat(z(ind, :), Nc*N, 1)).^2, 2) < eps_nn(i)^2);
    length(idx)
    
    plot3(zclouds2(idx, 1), zclouds2(idx, 2), zclouds2(idx, 3), '.', 'color', colors(i))
    hold on
    
    proj_data = zclouds2(idx, :) * proj_vec;
    
    C = cov(proj_data);
    
    [V, D] = eig(C);
    
    D_all(i, :) = sort(diag(D)', 'descend');
    
end

%%
for i=2:nvecs
D_all(:, i) = D_all(:, i) ./ D_all(:, 1);
end
figure; plot(D_all(:, 2:nvecs))



    
    

