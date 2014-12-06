clear all
close all

rng(321);

markersize = 500;

nsamples = 20;

D = 1;
lambda = linspace(1, 1000, nsamples);
s = sqrt(lambda / D);
neigs = 20;

mds_var = [];

for i=1:nsamples
    [V_tmp, D_tmp] = simulate_and_calculate_embedding(lambda(i), s(i));
    
    diff_dist = squareform(pdist(V_tmp * D_tmp));
    
    [~, E_tmp] = cmdscale(diff_dist);
    
    mds_var = [mds_var E_tmp];
    
end

%%
figure;
plot(lambda, mds_var(1,:))
hold on
plot(lambda, mds_var(2,:))
plot(lambda, mds_var(3,:))





