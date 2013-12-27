clear all
close all

fmt = '-djpeg';
res = '-r600';
set(0,'DefaultAxesFontSize',6)
set(0,'DefaultFigurePaperUnits','inches')

%% make figures for rxn example

species = {'S','E','E:S','S:E','D"S*','F:E*'};

% load data
load('../../rxn_example/snippet_data.mat','n', 'data');

% define data
N = size(data, 1);
train_ind = 1:1500;
test_ind = 1000:N;
[test_ind2, idx2] = setdiff(test_ind, train_ind);
[~, train_overlap_ind, test_overlap_ind] = intersect(train_ind, test_ind);

partial_ind = [1,2,4,6];
partial_ind2 = [1,2,3,5];

D = 2;
eps = 1e7;

% plot original data
az = [-5+180,-25+180,85,170];
el = [15,10,30,10];
az2 = [-10,-15,130,150];
el2 = [-55,-55,-35,-40];

for i=3:6
    figure;
    set(gcf,'paperposition',[0 0 3 2.5])
    scatter3(data(:,1),data(:,2),data(:,i),100,data(:,i),'.')
    xlabel(species{1})
    ylabel(species{2})
    zlabel(species{i})
    grid on
    view(az(i-2),el(i-2))
    saveas(gcf,sprintf('rxn_manifold%d',i-2),'fig')
    print(sprintf('rxn_manifold%d',i-2), fmt,res)
    
    set(gcf,'paperposition',[0 0 2 2])
    view(az2(i-2),el2(i-2))
    saveas(gcf,sprintf('rxn_manifold%d_2',i-2),'fig')
    print(sprintf('rxn_manifold%d_2',i-2), fmt,res)
end

% compute embeddings
inv_c = calc_covar(n(:,partial_ind,train_ind), D);
psi_mat = NL_ICA(data(train_ind,partial_ind), inv_c, eps);
%psi_mat = DM(data(train_ind,partial_ind), eps);

inv_c2 = calc_covar(n(:,partial_ind2,test_ind), D);
psi_mat2 = NL_ICA(data(test_ind,partial_ind2), inv_c2, eps);
%psi_mat2 = DM(data(test_ind,partial_ind2), eps);

for i=1:size(psi_mat,2)
    psi_mat(:,i) = sign(mean(psi_mat(:,i))) * psi_mat(:,i) / norm(psi_mat(:,i)) * size(psi_mat,1);
    psi_mat2(:,i) = sign(mean(psi_mat2(:,i))) * psi_mat2(:,i) / norm(psi_mat2(:,i)) * size(psi_mat2,1);
end

% plot embeddings
psi_coords = [1,3];

figure;
set(gcf,'paperposition',[0 0 3 2])
scatter(psi_mat(:,psi_coords(1)),psi_mat(:,psi_coords(2)),50,data(train_ind,1),'.')
xlabel(sprintf('\\psi_%d (data set 1)',psi_coords(1)))
ylabel(sprintf('\\psi_%d (data set 1)',psi_coords(2)))
axis([-100 250 -100 250])
saveas(gcf,'rxn_NLICA1','fig')
print('rxn_NLICA1',fmt,res)

figure;
set(gcf,'paperposition',[0 0 3 2])
scatter(psi_mat2(:,psi_coords(1)),psi_mat2(:,psi_coords(2)),50,data(test_ind,1),'.')
xlabel(sprintf('\\psi_%d (data set 2)',psi_coords(1)))
ylabel(sprintf('\\psi_%d (data set 2)',psi_coords(2)))
axis([-100 250 -100 250])
saveas(gcf,'rxn_NLICA2','fig')
print('rxn_NLICA2',fmt,res)

% plot correlation between embeddings

for i=1:length(psi_coords)
    figure;
    set(gcf,'paperposition',[0 0 3 2])
    plot(psi_mat(train_overlap_ind,psi_coords(i)), psi_mat2(test_overlap_ind, psi_coords(i)), '.')
    xlabel(sprintf('\\psi_%d (data set 1)',psi_coords(i)))
    ylabel(sprintf('\\psi_%d (data set 2)',psi_coords(i)))
    saveas(gcf,sprintf('rxn_NLICA_corr%d',i),'fig')
    print(sprintf('rxn_NLICA_corr%d',i),fmt,res)
end

corr(psi_mat(train_overlap_ind,psi_coords),psi_mat2(test_overlap_ind, psi_coords))

%% plot LapPyr
load Test_EstC4
figure;
set(gcf,'paperposition',[0 0 3 2])
plot(data(test_ind,4),Test_EstC4{end},'.')
xlabel(sprintf('true %s',species{4}))
ylabel(sprintf('reconstructed %s',species{4}))
saveas(gcf,'rxn_recon4','fig')
print('rxn_recon4',fmt,res)

load Test_EstC6
figure;
set(gcf,'paperposition',[0 0 3 2])
plot(data(test_ind,6),Test_EstC6{end},'.')
xlabel(sprintf('true %s',species{6}))
ylabel(sprintf('reconstructed %s',species{6}))
saveas(gcf,'rxn_recon6','fig')
print('rxn_recon6',fmt,res)


fprintf('nMSE for component 4: %2.4d \n',mean((data(test_ind,4)-Test_EstC4{end}).^2)/mean(data(test_ind,4).^2))
fprintf('nMSE for component 4: %2.4d \n',mean((data(test_ind,6)-Test_EstC6{end}).^2)/mean(data(test_ind,6).^2))

%% NN
K = 3;
knn_ind = knnsearch(data(1:999,1:2), data(1000:end,1:2), 'K', K);

Test_EstC4_nn = zeros(2001,1);
Test_EstC6_nn = zeros(2001,1);
for i=1:2001
    Test_EstC4_nn(i) = mean(data(knn_ind(i,:),4));
    Test_EstC6_nn(i) = mean(data(knn_ind(i,:),6));
end

figure;
plot(data(test_ind,4),Test_EstC4{end},'.')
hold on
plot(data(test_ind,4),Test_EstC4_nn,'.r')
xlabel(sprintf('true %s',species{4}))
ylabel(sprintf('reconstructed %s',species{4}))
legend('LapPyr','nn','location','best')
%saveas(gcf,'rxn_recon4','fig')
print('rxn_nn_compare4',fmt,res)

figure;
plot(data(test_ind,6),Test_EstC6{end},'.')
hold on
plot(data(test_ind,6),Test_EstC6_nn,'.r')
xlabel(sprintf('true %s',species{6}))
ylabel(sprintf('reconstructed %s',species{6}))
legend('LapPyr','nn','location','best')
%saveas(gcf,'rxn_recon6','fig')
print('rxn_nn_compare6',fmt,res)