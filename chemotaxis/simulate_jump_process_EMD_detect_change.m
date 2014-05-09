clear all
close all

rng(321);

markersize = 500;

%nsamples = 25;
nsamples = 20;
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

save('simulate_jump_process_EMD_detect_change_data')

%%
figure;
set(gcf, 'papersize', [8.6 4.3], 'units','centimeters')
set(gcf, 'paperposition', [0 0 8.6 4.3], 'units','centimeters')
errorbar(lambda, mean(D2), std(D2), '.b', 'markersize', 14)
hold on
errorbar(lambda, mean(D3), std(D3), '.r', 'markersize', 14)
xlabel('$\lambda$', 'interpreter','latex', 'fontsize', 20)
ylabel('$\mu_k$', 'interpreter','latex', 'fontsize', 20)
axis([-20 420 0.37 0.57])
legend('k = 1','k = 2','location','northwest')

saveas(gcf, 'detect_change_eigenvalues', 'epsc')



