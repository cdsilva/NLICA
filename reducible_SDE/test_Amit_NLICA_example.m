% this code simulates the data from "Non LInear INdependent Component
% Analysis"

clear all
close all

res = '-r300';
fmt = '-djpeg';

rng(321);

%% define parameters

% number of data points
N = 2000;

% dimension
dim = 2;

% number of points in each cloud
Nc = 200;

% dt for clouds
dt1 = 0.001;

%% do simulations

x = rand(N, dim);

f = @(x) [x(:,1)+x(:,2).^3 x(:,2)-x(:,1).^3];
y = f(x);

%% plot results

figure;
plot(x(:,1),x(:,2),'.')
xlabel('x_1')
ylabel('x_2')
print('xdata', fmt, res)

figure;
plot(y(:,1),y(:,2),'.')
xlabel('y_1')
ylabel('y_2')
print('ydata', fmt, res)

%% calculate clouds

xclouds = dt1 * randn(Nc, dim, N);
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
xlabel('x_1')
ylabel('x_2')
print('xdata_colored_NIV1', fmt, res)

figure;
scatter(x(:,1),x(:,2),200,V(:,3),'.')
xlabel('x_1')
ylabel('x_2')
print('xdata_colored_NIV2', fmt, res)

%% dmaps
W = squareform(pdist(y)).^2;
[V, D] = dmaps(W, median(W(:)), neigs);

figure;
scatter(x(:,1),x(:,2),200,V(:,2),'.')
print('xdata_colored_DMAPS1', fmt, res)

figure;
scatter(x(:,1),x(:,2),200,V(:,3),'.')
print('xdata_colored_DMAPS2', fmt, res)

%% do simulation with noise
dim = 3;
dt2 = dt1 / 100;

z = [y zeros(N, 1)];

zclouds = dt2 * randn(Nc, dim, N);
for i=1:N
    zclouds(:,1:2,i) = zclouds(:,1:2,i) + yclouds(:,:,i);
end


% calculate inverse covariances
[inv_c, ~, ranks] = covariances2(zclouds, dim);

% NIV
if size(z,1) > 4000
    disp('Too much data; NIV will fail')
    return
end
[V, D] = NIV(z, inv_c, 0, neigs, 0);

figure;
scatter(x(:,1),x(:,2),200,V(:,2),'.')
xlabel('x_1')
ylabel('x_2')
print('xdata_noise1_colored_NIV1', fmt, res)

figure;
scatter(x(:,1),x(:,2),200,V(:,3),'.')
xlabel('x_1')
ylabel('x_2')
print('xdata_noise1_colored_NIV2', fmt, res)

%% do simulation with noise

dt2 = dt1 * 100;

z = [y zeros(N, 1)];

zclouds = dt2 * randn(Nc, dim, N);
for i=1:N
    zclouds(:,1:2,i) = zclouds(:,1:2,i) + yclouds(:,:,i);
end

% calculate inverse covariances
[inv_c, ~, ranks] = covariances2(zclouds, dim);

% NIV
if size(z,1) > 4000
    disp('Too much data; NIV will fail')
    return
end
[V, D] = NIV(z, inv_c, 0, neigs, 0);

figure;
scatter(x(:,1),x(:,2),200,V(:,2),'.')
xlabel('x_1')
ylabel('x_2')
print('xdata_noise2_colored_NIV1', fmt, res)

figure;
scatter(x(:,1),x(:,2),200,V(:,3),'.')
xlabel('x_1')
ylabel('x_2')
print('xdata_noise2_colored_NIV2', fmt, res)
