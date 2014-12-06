clear all
close all

rng(123);


%%

ntraj = 10;
npoints = 1000;
x = repmat(linspace(0, 1, npoints), ntraj, 1);

y = zeros(size(x));
for i=1:ntraj
    y(i, :) = randn * exp(-10*x(i,:));
end

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
y = randn(size(x)) .* (1+8*exp(-4*x));

make_fig(3, 3);
scatter(x, smooth(y, 11), 50, -x,'.')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca, 'ylim', [-8 8])
% axis off
saveas(gcf, 'schematic_DS2.eps', 'epsc');


