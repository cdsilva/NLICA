clear all
close all

%%
% 
% N = 1500; % number of points considered
% t = rand(1,N);
% t = sort(4*pi*sqrt(t))'; 
% 
% h = 15*pi;
% 
% z = h*rand(N,1); % random heights
% x = (t+.1).*cos(t);
% y = (t+.1).*sin(t);
% data = [x,y,z]; % data
% 
% %%
% figure;
% plot3(x,y,z,'.')
% 
% axis equal
% 
% 
% %%
% 
% W = squareform(pdist(data)).^2;
% eps = 5;
% 
% neigs = 10;
% 
% [V, D] = dmaps(W, eps, neigs);
% 
% 
% %%
% figure;
% scatter3(x,y,z,50,V(:,2),'.')
% 
% figure;
% scatter3(x,y,z,50,V(:,3),'.')
% 
% %%
% 
% eps_med_scale = 3;
% res = compute_residuals_DMAPS(V, eps_med_scale);
% 
% figure;
% plot(res, 'o-')

%%

load swissroll1.mat

make_fig(2,2);
plot3(x,y,z,'.', 'markersize', 10)
view(-20, 70)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
print('swissroll1.eps','-depsc')

% 
% figure;
% bar(res(2:end))
% ylabel('cross-validation error')
% print('swissroll1_cv.eps','-depsc')
% 
make_fig(3,2);
scatter3(x,y,z,50,sign(corr(V(:,2),t))*V(:,2),'.')
view(-20, 75)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_1');
print('swissroll1_color1.eps','-depsc')

make_fig(3,2);
scatter3(x,y,z,50,sign(corr(V(:,3),z))*V(:,3),'.')
view(-20, 75)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_2');
print('swissroll1_color2.eps','-depsc')

make_fig(3,2);
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'ylim', [0.9 1])
xlabel('k')
ylabel('\mu_k')
axis square
print('swissroll1_evals.eps','-depsc')

%%

figure;

D1 = pdist(V(:,2:3));
D2 = pdist(V(:,2:4));

subplot(2,2,1)
plot(V(:,2),V(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')

eps = logspace(-8, -3, 100);
sumA1 = zeros(size(eps));
sumA2 = zeros(size(eps));

for i=1:length(eps)
    sumA1(i) = sum(exp(-D1.^2/(2*eps(i))));
    sumA2(i) = sum(exp(-D2.^2/(2*eps(i))));
end

subplot(2,2,2)
loglog(eps, sumA1, '.')
hold on
loglog(eps, 1e10*eps)
p = polyfit(log(eps), log(sumA1), 1);
legend(sprintf('data; est dim = %2.2f', p(1)*2), 'dim=2 line','location','best')
xlabel('\epsilon')
ylabel('$\sum W_{ij}$', 'interpreter','latex')

subplot(2,2,3)
plot3(V(:,2),V(:,3),V(:,4),'.')
xlabel('\phi_2')
ylabel('\phi_3')
zlabel('\phi_4')

subplot(2,2,4)
loglog(eps, sumA2, '.r')
hold on
loglog(eps, 1e10*eps, 'r')
p = polyfit(log(eps), log(sumA2), 1);
legend(sprintf('data; est dim = %2.2f', p(1)*2), 'dim=2 line','location','best')
xlabel('\epsilon')
ylabel('$\sum W_{ij}$', 'interpreter','latex')

%%

load swissroll2.mat

make_fig(2,2);
plot3(x,y,z,'.', 'markersize', 10)
view(-20, 75)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
print('swissroll2.eps','-depsc')
% 
% figure;
% bar(res(2:end))
% ylabel('cross-validation error')
% print('swissroll2_cv.eps','-depsc')
% 

make_fig(3,2);
scatter3(x,y,z,50,sign(corr(V(:,2),t))*V(:,2),'.')
view(-20, 75)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_1');
print('swissroll2_color1.eps','-depsc')

make_fig(3,2);
scatter3(x,y,z,50,sign(corr(V(:,6),z))*V(:,6),'.')
view(-20, 75)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_5');
print('swissroll2_color2.eps','-depsc')


make_fig(3,2);
make_colored_bars(diag(D(2:end,2:end)), res(2:end))
set(gca, 'ylim', [0.9 1])
xlabel('k')
ylabel('\mu_k')
axis square
print('swissroll2_evals.eps','-depsc')

%%

figure;

D1 = pdist(V(:,2:3));
D2 = pdist(V(:,[2 6]));

subplot(2,2,1)
plot(V(:,2),V(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')

eps = logspace(-8, -4, 100);
sumA1 = zeros(size(eps));
sumA2 = zeros(size(eps));

for i=1:length(eps)
    sumA1(i) = sum(exp(-D1.^2/(2*eps(i))));
    sumA2(i) = sum(exp(-D2.^2/(2*eps(i))));
end

subplot(2,2,2)
loglog(eps, sumA1, '.')
hold on
loglog(eps, 1e7*eps.^(0.5))
p = polyfit(log(eps), log(sumA1), 1);
legend(sprintf('data; est dim = %2.2f', p(1)*2), 'dim=1 line','location','best')
xlabel('\epsilon')
ylabel('$\sum W_{ij}$', 'interpreter','latex')

subplot(2,2,3)
plot(V(:,2),V(:,6),'.')
xlabel('\phi_2')
ylabel('\phi_6')

subplot(2,2,4)
loglog(eps, sumA2, '.r')
hold on
loglog(eps, 1e10*eps, 'r')
p = polyfit(log(eps), log(sumA2), 1);
legend(sprintf('data; est dim = %2.2f', p(1)*2), 'dim=2 line','location','best')
xlabel('\epsilon')
ylabel('$\sum W_{ij}$', 'interpreter','latex')




