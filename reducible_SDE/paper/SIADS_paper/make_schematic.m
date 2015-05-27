clear all
close all

rng(123);


%%

ntraj = 2;
npoints = 1000;
t = repmat(linspace(0, 1, npoints), ntraj, 1);

x1 = exp(-2*t);
x2 = exp(-50*t);


x1(1,:) = -0.95*x1(1,:);
x1(2,:) = -0.7*x1(2,:);

x2(1,:) = 0.75*x2(1,:);
x2(2,:) = -0.65*x2(2,:);


% figure;
% plot(x, y, '.')
% set(gca, 'xtick', [])
% set(gca, 'ytick', [])

make_fig(3, 3);
scatter(x1(:), x2(:), 20, t(:),'.')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca, 'ylim', [-1 1])
xlabel('x_1')
ylabel('x_2')
% axis off
colormap(flipud(jet))
h = colorbar;
xlabel(h,'t');
set(h,'YTick',[0 1])
saveas(gcf, 'schematic_DS1.eps', 'epsc');

% figure;
% plot(-(y'+1).*cos((pi/2)*x'), (y'+1).*sin((pi/2)*x'))

%%

npoints = 4000;

x = linspace(0, 1, npoints);
y = randn(size(x)) .* max(10*exp(-8*x), 2);

make_fig(3, 3);
scatter(x, smooth(y, 11), 50, x/max(x),'.')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca, 'ylim', [-3 3])
xlabel('x_1')
ylabel('x_2')
% axis off
colormap(flipud(jet))
h = colorbar;
xlabel(h,'t');
set(h,'YTick',[0 1])
saveas(gcf, 'schematic_DS2.eps', 'epsc');


