clear all
close all

rng(321);

markersize = 500;

%nsamples = 25;
nsamples = 5;
nreplicates = 5;

D = 1;
lambda = linspace(1, 2500, nsamples);
s = sqrt(lambda / D);
neigs = 4;


D_all = zeros(nreplicates, nsamples, neigs);
D2_all = zeros(nreplicates, nsamples, neigs);
res = zeros(nreplicates, nsamples, neigs);
for i=1:nsamples
    for j=1:nreplicates
        [V_tmp, D_tmp] = simulate_and_calculate_embedding(lambda(i), s(i));
        D_all(j,i,:) = diag(D_tmp(1:neigs,1:neigs));
        [V2, D2] = reorder_DMAPS_noharmonics(V_tmp(:,1:neigs), D_tmp(1:neigs,1:neigs));
        D2_all(j,i,:) = diag(D2);
        res(j,i,:) = compute_residuals_DMAPS(V_tmp(:,1:neigs), D_tmp(1:neigs,1:neigs));

    end
end

%%
figure;
errorbar(lambda, mean(D_all(:,:,2)), std(D_all(:,:,2)))
hold on
errorbar(lambda, mean(D_all(:,:,3)), std(D_all(:,:,3)))
errorbar(lambda, mean(D_all(:,:,4)), std(D_all(:,:,4)))
xlabel('\lambda')
ylabel('eigenvalues')

figure;
errorbar(lambda, mean(D2_all(:,:,2)), std(D2_all(:,:,2)))
hold on
errorbar(lambda, mean(D2_all(:,:,3)), std(D2_all(:,:,3)))
errorbar(lambda, mean(D2_all(:,:,4)), std(D2_all(:,:,4)))
xlabel('\lambda')
ylabel('adjusted eigenvalues')

figure;
plot(lambda, mean(D_all(:,:,2)) - mean(D_all(:,:,3)))

figure;
plot(lambda, mean(D2_all(:,:,2)) - mean(D2_all(:,:,3)))

% save('simulate_jump_process_EMD_detect_change_data')
%
% %%
% figure;
% set(gcf, 'papersize', [8.6 4.3], 'units','centimeters')
% set(gcf, 'paperposition', [0 0 8.6 4.3], 'units','centimeters')
% errorbar(lambda, mean(D2), std(D2), '.b', 'markersize', 14)
% hold on
% errorbar(lambda, mean(D3), std(D3), '.r', 'markersize', 14)
% xlabel('$\lambda$', 'interpreter','latex', 'fontsize', 20)
% ylabel('$\mu_k$', 'interpreter','latex', 'fontsize', 20)
% axis([-20 420 0.37 0.57])
% legend('k = 1','k = 2','location','northwest')
%
% saveas(gcf, 'detect_change_eigenvalues', 'epsc')

%%

D_rescaled
res_thres = 0.4;

for i=1:nsamples
    for j=1:nreplicates
        [V_tmp, D_tmp] = simulate_and_calculate_embedding(lambda(i), s(i));
        D_all(j,i,:) = diag(D_tmp(1:neigs,1:neigs));
        [V2, D2] = reorder_DMAPS_noharmonics(V_tmp(:,1:neigs), D_tmp(1:neigs,1:neigs));
        D2_all(j,i,:) = diag(D2);
        res(j,i,:) = compute_residuals_DMAPS(V_tmp(:,1:neigs), D_tmp(1:neigs,1:neigs));

    end
end




