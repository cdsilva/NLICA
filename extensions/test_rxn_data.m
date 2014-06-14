clear all
close all

load ../rxn_example/for_lp.mat

data = orig_data;

figure;
plot3(data(:,1),data(:,2),data(:,3),'.')

%%
% idx = find(data(:,1) > 5e4);
subregion_idx = find(psi_mat_full(:,1).^2 + psi_mat_full(:,3).^2 < 0.02^2);
data2 = data(subregion_idx, :);
figure;
plot3(data(:,1),data(:,2),data(:,3),'.')
hold on
plot3(data2(:,1),data2(:,2),data2(:,3),'.r')
[~, idx_edge] = min(data2(:, 1));
plot3(data2(idx_edge,1),data2(idx_edge,2),data2(idx_edge,3),'.g', 'markersize', 20)

%%
K = 50;
IDX = knnsearch(data2,data2(idx_edge,:), 'K', K);

figure;
plot3(data(:,1),data(:,2),data(:,3),'.')
hold on
plot3(data2(:,1),data2(:,2),data2(:,3),'.r')
plot3(data2(IDX,1),data2(IDX,2),data2(IDX,3),'.g', 'markersize', 20)


%%

delta = data2(idx_edge,:) - mean(data2(IDX,:));

new_point = data2(idx_edge,:) + delta;

figure;
plot3(data(:,1),data(:,2),data(:,3),'.')
hold on
plot3(data2(:,1),data2(:,2),data2(:,3),'.r')
plot3(new_point(1), new_point(2), new_point(3),'.g', 'markersize', 20)

%%

figure; 
plot(psi_mat_full(:,1),psi_mat_full(:,3),'.')
hold on
plot(psi_mat_full(subregion_idx,1),psi_mat_full(subregion_idx,3),'.r')

%%

embed_coord = psi_mat_full(:, [1 3]);

save('data_for_elio.mat', 'data', 'embed_coord', 'subregion_idx');