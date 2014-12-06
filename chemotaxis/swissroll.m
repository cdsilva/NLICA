clear all
close all

%%

N = 1500; % number of points considered
t = rand(1,N);
t = sort(4*pi*sqrt(t))'; 

h = 15*pi;

z = h*rand(N,1); % random heights
x = (t+.1).*cos(t);
y = (t+.1).*sin(t);
data = [x,y,z]; % data

%%
figure;
plot3(x,y,z,'.')

axis equal


%%

W = squareform(pdist(data)).^2;
eps = 5;

neigs = 10;

[V, D] = dmaps(W, eps, neigs);


%%
figure;
scatter3(x,y,z,50,V(:,2),'.')

figure;
scatter3(x,y,z,50,V(:,3),'.')

%%

eps_med_scale = 3;
res = compute_residuals_DMAPS(V, eps_med_scale);

figure;
plot(res, 'o-')

%%

load swissroll1.mat

figure;
plot3(x,y,z,'.')
axis equal
print('swissroll1.eps','-depsc')

figure;
bar(res(2:end))
ylabel('cross-validation error')
print('swissroll1_cv.eps','-depsc')

figure;
scatter3(x,y,z,50,V(:,2),'.')
axis equal
print('swissroll1_color1.eps','-depsc')

figure;
scatter3(x,y,z,50,V(:,3),'.')
axis equal
print('swissroll1_color2.eps','-depsc')

%%

load swissroll2.mat

figure;
plot3(x,y,z,'.')
axis equal
print('swissroll2.eps','-depsc')

figure;
bar(res(2:end))
ylabel('cross-validation error')
print('swissroll2_cv.eps','-depsc')

figure;
scatter3(x,y,z,50,V(:,2),'.')
axis equal
print('swissroll2_color1.eps','-depsc')

figure;
scatter3(x,y,z,50,V(:,6),'.')
axis equal
print('swissroll2_color2.eps','-depsc')





