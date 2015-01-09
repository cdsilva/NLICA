clear all
close all

rng(123);


%%

ntraj = 2;
npoints = 1000;
x = repmat(linspace(0, 1, npoints), ntraj, 1);

y = zeros(size(x));
y(1, :) = exp(-10*x(1,:));
y(2, :) = -exp(-20*x(2,:));

% figure;
% plot(x, y, '.')
% set(gca, 'xtick', [])
% set(gca, 'ytick', [])

make_fig(3, 3);
scatter(x(:), y(:), 50, x(:)/max(x(:)),'.')
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


