clear all
close all

rng(321);

markersize = 500;

%nsamples = 25;
nsamples = 10;
nreplicates = 10;

D = 1;
lambda = linspace(1, 400, nsamples);
s = sqrt(lambda / D);

D2 = zeros(nreplicates, nsamples);
D3 = zeros(nreplicates, nsamples);
for i=1:nsamples
    for j=1:nreplicates
        [V2, D_tmp] = simulate_and_calculate_embedding(lambda(i), s(i));
        D2(j, i) = D_tmp(2,2);
        D3(j, i) = D_tmp(3,3);
    end
end

%%
figure;
errorbar(lambda, mean(D2), std(D2))
hold on
errorbar(lambda, mean(D3), std(D3))

% scatter(lambda, D2(2,:), markersize, 'b', '.')
% hold on
% scatter(lambda, D2(3,:), markersize, 'r', '.')
% xlabel('\lambda')
% ylabel('\mu_k')
% legend('k=1','k=2','location','best')
% 

