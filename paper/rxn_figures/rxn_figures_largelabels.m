clear all
close all

fmt = '-djpeg';
res = '-r600';
set(0,'DefaultAxesFontSize',6)
set(0,'DefaultFigurePaperUnits','inches')

%% make figures for rxn example

species = {'S','E','E:S','S:E','D:S*','F:E*'};

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
    xlabel(species{1},'fontsize',14)
    ylabel(species{2},'fontsize',14)
    zlabel(species{i},'fontsize',14)
    grid on
    view(az(i-2),el(i-2))
    %saveas(gcf,sprintf('rxn_manifold%d',i-2),'fig')
    print(sprintf('rxn_manifold_largelabels_%d',i-2), fmt,res)
    
%     set(gcf,'paperposition',[0 0 2 2])
%     view(az2(i-2),el2(i-2))
%     saveas(gcf,sprintf('rxn_manifold%d_2',i-2),'fig')
%     print(sprintf('rxn_manifold%d_2',i-2), fmt,res)
end


%% LapPyr
load Test_EstC4
figure;
set(gcf,'paperposition',[0 0 3 2])
plot(data(test_ind,4),Test_EstC4{end},'.')
xlabel(sprintf('true %s',species{4}),'fontsize',14)
ylabel(sprintf('reconstructed %s',species{4}),'fontsize',14)
print('rxn_recon4_largelabels',fmt,res)

load Test_EstC6
figure;
set(gcf,'paperposition',[0 0 3 2])
plot(data(test_ind,6),Test_EstC6{end},'.')
xlabel(sprintf('true %s',species{6}),'fontsize',14)
ylabel(sprintf('reconstructed %s',species{6}),'fontsize',14)
print('rxn_recon6_largelabels',fmt,res)