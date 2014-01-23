clear all
close all

dim = 2;
n = 1000;
scale = [1 0.1];

data = rand(n, dim);
for i=1:dim
    data(:,i) = data(:,i) * scale(i);
end

figure; 
plot(data(:,1),data(:,2),'.')
axis equal
xlabel('x')
ylabel('y')

W = squareform(pdist(data)).^2;
eps = median(W(:));

neigs = 50;
[V, D] = dmaps(W, eps, neigs);

dcorr_cum = zeros(neigs-2, 1);

for i=3:neigs
    dcorr_cum(i-2) = dist_corr(V(:,2:i-1), V(:,i));
end

figure; 
plot(3:neigs, dcorr_cum, '.', 'markersize', 12)
xlabel('k')
ylabel('distance correlation between \phi_k and \phi_2-\phi_{k-1}')

[~, idx] = min(dcorr_cum);
idx = idx + 2;

figure;
plot(V(:,2),V(:,idx),'.')
xlabel('\phi_2')
ylabel(sprintf('\\phi_%d', idx))
