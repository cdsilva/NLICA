clear all
close all

fmt = '-depsc';
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
test_ind = 1:2000;

partial_ind = [1,2,4,6];

D = 2;
eps = 1e7;

npoints = [10 100 250 500 1000 1500];
[inv_c, spec] = calc_covar2(n(:,partial_ind,1:npoints(end)), D);
%f1 = figure;
%set(gcf,'paperposition',[0 0 5 8])
%f2 = figure;
%set(gcf,'paperposition',[0 0 5 8])
for j=1:length(npoints)
    train_ind = 1:npoints(j);
    
    % compute embeddings
    [psi_mat, D] = NL_ICA2(data(train_ind,partial_ind), inv_c(:,:,train_ind), eps);
    
    for i=1:size(psi_mat,2)
        psi_mat(:,i) = sign(mean(psi_mat(:,i))) * psi_mat(:,i) / norm(psi_mat(:,i)) * size(psi_mat,1);
    end
    
    if psi_mat(1,1) > 0
        psi_mat(:,1) = -psi_mat(:,1);
    end
    if psi_mat(1,3) < 0
        psi_mat(:,3) = -psi_mat(:,3);
    end
    
    % plot embeddings
    psi_coords = [1,3];
    
    figure;
    set(gcf,'paperposition',[0 0 1.5 1.25])
    scatter(psi_mat(:,psi_coords(1)),psi_mat(:,psi_coords(2)),50,data(train_ind,1),'.')
    xlabel(sprintf('\\psi_%d',psi_coords(1)))
    ylabel(sprintf('\\psi_%d',psi_coords(2)))
    %title(sprintf('n = %d',npoints(j)))
    saveas(gcf,sprintf('rxn_npoints_embed%d',j),'fig')
    print(sprintf('rxn_npoints_embed%d',j), fmt,res)
    
    figure;
    set(gcf,'paperposition',[0 0 1.5 1.25])
    plot(diag(D(2:end,2:end)),'.','markersize',10)
    ylim([0 0.5e-3])
    xlabel('k')
    ylabel('\lambda_k')
    %title(sprintf('n = %d',npoints(j)))
    saveas(gcf,sprintf('rxn_npoints_spectrum%d',j),'fig')
    print(sprintf('rxn_npoints_spectrum%d',j), fmt,res)
end
% figure(f1)
% print('rxn_compare_npoints',fmt,res)
% figure(f2)
% print('rxn_compare_npoints2',fmt,res)

