clear all 
close all

fmt = '-djpeg';
res = '-r600';
set(0,'DefaultAxesFontSize',6)
set(0,'DefaultFigurePaperUnits','inches')

%% Load Data
%load '~/resesarch/Alanine/Short_Trajs_Cartesian_Coord.dat'
%load '~/resesarch/Alanine/Short_Trajs_dihedrals.dat'
load '../../Ronen/Short_Trajs_Cartesian_Coord.dat'
load '../../Ronen/Short_Trajs_dihedrals.dat'

time_indices1 = 50001:60000;
time_indices2 = 58001:70000;
num_overlap = 2000;

coor_indices1 = 2:31; % all
coor_indices2 = [2,3,4,8,9,10,14,15,16,20,21,22,26,27,28]; % only odd

%% Odd Coordinates
data = Short_Trajs_Cartesian_Coord(time_indices1, coor_indices1);
    
% Compute local covariances
covariances;
    
% NLICA
NLICA;
psi_mat_nlica1 = psi_mat;

% DM
DM;   
psi_mat_dm1 = psi_mat;

%% Even Coordinates  
data = Short_Trajs_Cartesian_Coord(time_indices2, coor_indices2);
    
% Compute local covariances
covariances;
    
% NLICA
NLICA;
psi_mat_nlica2 = psi_mat;

% DM
DM;   
psi_mat_dm2 = psi_mat;

%% align two embeddings
for i=1:3
    psi_mat_nlica1(:,i) = psi_mat_nlica1(:,i) / norm(psi_mat_nlica1(:,i)) * size(psi_mat_nlica1,1);
    psi_mat_nlica2(:,i) = psi_mat_nlica2(:,i) / norm(psi_mat_nlica2(:,i)) * size(psi_mat_nlica2,1);
    corr1 = sum(psi_mat_nlica1(end-num_overlap+1:end,i).*psi_mat_nlica2(1:num_overlap,i)) / sum(psi_mat_nlica1(end-num_overlap+1:end,i).^2);
    if corr1 < 0
        psi_mat_nlica2(:,i) = -psi_mat_nlica2(:,i);
    end
    psi_mat_dm1(:,i) = psi_mat_dm1(:,i) / norm(psi_mat_dm1(:,i)) * size(psi_mat_dm1,1);
    psi_mat_dm2(:,i) = psi_mat_dm2(:,i) / norm(psi_mat_dm2(:,i)) * size(psi_mat_dm2,1);
    corr1 = sum(psi_mat_dm1(end-num_overlap+1:end,i).*psi_mat_dm2(1:num_overlap,i)) / sum(psi_mat_dm1(end-num_overlap+1:end,i).^2);
    if corr1 < 0
        psi_mat_dm2(:,i) = -psi_mat_dm2(:,i);
    end
end

disp('DM correlation:')
corr(psi_mat_dm1(end-num_overlap+1:end,1:3),psi_mat_dm2(1:num_overlap,1:3))

disp('NLICA correlation:')
corr(psi_mat_nlica1(end-num_overlap+1:end,1:3),psi_mat_nlica2(1:num_overlap,1:3))

% plot correlations
for i=1:3
    figure;
    plot(psi_mat_nlica1(end-num_overlap+1:end,i),psi_mat_nlica2(1:num_overlap,i),'.')
    xlabel(sprintf('\\psi_%d^{all}',i))
    ylabel(sprintf('\\psi_%d^{odd}',i))
    %print(sprintf('ala2_corr%d',i),fmt,res)
end


for i=1:3
    figure;
    plot(psi_mat_dm1(end-num_overlap+1:end,i),psi_mat_dm2(1:num_overlap,i),'.')
    xlabel(sprintf('DM%d^{all}',i))
    ylabel(sprintf('DM%d^{odd}',i))
    %print(sprintf('ala2_corr%d',i),fmt,res)
end

%% 3D embeddings

arrow_start = -250;
arrow_len = 100;
arrows_color = eye(3);

order = [1 2 3];
figure; 
set(gcf,'paperposition',[0 0 2 2])
scatter3(psi_mat_nlica1(:,order(1)),psi_mat_nlica1(:,order(2)),psi_mat_nlica1(:,order(3)),10, Short_Trajs_Cartesian_Coord(time_indices1, 3),'.')
hold on
for i=1:3
    length_arrow = zeros(3,1);
    length_arrow(i) = arrow_len;
    quiver3(arrow_start, arrow_start, arrow_start, length_arrow(1),length_arrow(2),length_arrow(3),'color', arrows_color(order(i),:))
end
grid on
xlabel(sprintf('\\psi_%d',order(1)))
ylabel(sprintf('\\psi_%d',order(2)))
zlabel(sprintf('\\psi_%d',order(3)))
axis([-inf inf -inf inf -inf inf])
view(-10,25)
saveas(gcf,'ala2_embed1','fig')
print('ala2_embed1',fmt,res)

order = [3 1 2];
figure; 
set(gcf,'paperposition',[0 0 2 2])
scatter3(psi_mat_nlica1(:,order(1)),psi_mat_nlica1(:,order(2)),psi_mat_nlica1(:,order(3)),10, Short_Trajs_dihedrals(time_indices1, 2),'.')
hold on
for i=1:3
    length_arrow = zeros(3,1);
    length_arrow(i) = arrow_len;
    quiver3(arrow_start, arrow_start, arrow_start, length_arrow(1),length_arrow(2),length_arrow(3), 'color', arrows_color(order(i),:))
end
grid on
xlabel(sprintf('\\psi_%d',order(1)))
ylabel(sprintf('\\psi_%d',order(2)))
zlabel(sprintf('\\psi_%d',order(3)))
axis([-inf inf -inf inf -inf inf])
view(-10,25)
saveas(gcf,'ala2_embed2','fig')
print('ala2_embed2',fmt,res)

order = [2 3 1];
figure; 
set(gcf,'paperposition',[0 0 2 2])
scatter3(psi_mat_nlica1(:,order(1)),psi_mat_nlica1(:,order(2)),psi_mat_nlica1(:,order(3)),10, Short_Trajs_dihedrals(time_indices1, 3),'.')
hold on
for i=1:3
    length_arrow = zeros(3,1);
    length_arrow(i) = arrow_len;
    quiver3(arrow_start, arrow_start, arrow_start, length_arrow(1),length_arrow(2),length_arrow(3), 'color', arrows_color(order(i),:))
end
grid on
xlabel(sprintf('\\psi_%d',order(1)))
ylabel(sprintf('\\psi_%d',order(2)))
zlabel(sprintf('\\psi_%d',order(3)))
view(-10,25)
axis([-inf inf -inf inf -inf inf])
saveas(gcf,'ala2_embed3','fig')
print('ala2_embed3',fmt,res)


%% Lap Pyr

load ala_data_forLP
load LP_A
load NN_A
test_ind = 1:5000;

% comps = 1:3;
% train_data = psi_mat_nlica1(:,comps);
% test_data = psi_mat_nlica2(num_overlap+1:end,comps);
% train_fn = Short_Trajs_Cartesian_Coord(time_indices1, 2:end);
% test_fn = Short_Trajs_Cartesian_Coord(time_indices2(num_overlap+1:end), 2:end);
% 
% save('ala_data_forLP.mat','train_data','test_data','train_fn','test_fn')
% 
% err_thresh = 1e-8;
% eps_lp = 10000;
% 
% %train_ind = randsample(size(train_data,1),2000);
% %test_ind = randsample(size(test_data,1),1000);
% train_ind = ceil(linspace(1,size(train_data,1),2000));
% test_ind = ceil(linspace(1,size(test_data,1),1000));
% 
% [train_est, test_est, ~, ~, ~, ~] = LP_chm(train_data(train_ind,:), test_data(test_ind,:), train_fn(train_ind,:), err_thresh, eps_lp);
% 
% true_pos = test_fn(test_ind,:);
% recon_pos = test_est{end};
% save('LapPyr_recon.mat','true_pos','recon_pos');

for i=1:30
	figure;
    set(gcf,'paperposition',[0 0 2 2])
    %plot(test_est{end}(:,i),test_fn(test_ind,i),'.')
    plot(test_fn(test_ind,i),LP_A(:,i),'.','markersize',1)
    avg1 = mean(test_fn(:,i));
    delta = 0.75;
    axis([avg1-delta avg1+delta avg1-delta avg1+delta])
    xlabel('true position')
    ylabel('reconstructed position')
    saveas(gcf,sprintf('ala2_recon_corr_%d_%d',ceil(i/3), mod(i-1,3)+1),'fig')
    print(sprintf('ala2_recon_corr_%d_%d',ceil(i/3), mod(i-1,3)+1),fmt,res)
end

bonds = [1 2;
    2 3; 
    2 4;
    4 5;
    5 6;
    5 7;
    7 8;
    7 9;
    9 10];
avg_bond = 0;
for i=1:9
    avg_bond = avg_bond + mean(sum((test_fn(test_ind,3*bonds(i,1)-2:3*bonds(i,1))-test_fn(test_ind,3*bonds(i,2)-2:3*bonds(i,2))).^2,2));
end
avg_bond = avg_bond / 9;

mse_atoms = zeros(10,1);
mse_atoms_nn = zeros(10,1);
%nmse_atoms = zeros(10,1);
%nmse_atoms_nn = zeros(10,1);
for i=1:3:30
    %mse_atoms(ceil(i/3)) = sum(sum((test_est{end}(:,i:i+2)-test_fn(test_ind,i:i+2)).^2))/length(test_ind);
    mse_atoms(ceil(i/3)) = sum(sum((LP_A(:,i:i+2)-test_fn(test_ind,i:i+2)).^2))/length(test_ind);
    mse_atoms_nn(ceil(i/3)) = sum(sum((NN_A(:,i:i+2)-test_fn(test_ind,i:i+2)).^2))/length(test_ind);
    %nmse_atoms(ceil(i/3)) = sum(sum((LP_A(:,i:i+2)-test_fn(test_ind,i:i+2)).^2))/sum(sum((test_fn(test_ind,i:i+2)-repmat(mean(test_fn(test_ind,i:i+2)),length(test_ind),1)).^2));
    %nmse_atoms_nn(ceil(i/3)) = sum(sum((NN_A(:,i:i+2)-test_fn(test_ind,i:i+2)).^2))/sum(sum((test_fn(test_ind,i:i+2)-repmat(mean(test_fn(test_ind,i:i+2)),length(test_ind),1)).^2));
end

figure;
set(gcf,'paperposition',[0 0 3 2])
bar([mse_atoms/avg_bond mse_atoms_nn/avg_bond])
xlabel('atom index')
ylabel('normalized mean squared reconstruction error')
legend('LapPyr','nn','location','northwest')
saveas(gcf,'ala2_recon_mse','fig')
print('ala2_recon_mse',fmt,res)

% figure;
% set(gcf,'paperposition',[0 0 3 2])
% bar([nmse_atoms nmse_atoms_nn])
% xlabel('atom index')
% ylabel('average reconstruction error')
% legend('LapPyr','nn','location','northwest')
%saveas(gcf,'ala2_recon_nmse','fig')
%print('ala2_recon_nmse',fmt,res)