clear all
close all

rng(123);


%%

ntraj = 2;
npoints = 1000;
x = repmat(linspace(0, 1, npoints), ntraj, 1);

y = zeros(size(x));
y(1, :) = exp(-20*x(1,:));
y(2, :) = -exp(-20*x(2,:));

% figure;
% plot(x, y, '.')
% set(gca, 'xtick', [])
% set(gca, 'ytick', [])

make_fig(3, 3);
scatter(x(:), y(:), 50, -x(:),'.')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca, 'ylim', [-1 1])
% axis off
saveas(gcf, 'schematic_DS1.eps', 'epsc');

% figure;
% plot(-(y'+1).*cos((pi/2)*x'), (y'+1).*sin((pi/2)*x'))

%%

npoints = 4000;

x = linspace(0, 1, npoints);
y = randn(size(x)) .* max(10*exp(-7*x), 2);

make_fig(3, 3);
scatter(x, smooth(y, 11), 50, -x,'.')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca, 'ylim', [-3 3])
% axis off
saveas(gcf, 'schematic_DS2.eps', 'epsc');


