clear all
close all

%%

fmt = '-djpeg';
res = '-r600';

load data_for_elio.mat

load EDGES.mat

data2 = data(subregion_idx,:);
embed_coord2 = embed_coord(subregion_idx,:);

%%

make_fig;
plot(embed_coord(:,1), embed_coord(:,2),'.', 'color', 0.9*ones(1,3))
hold on
plot(embed_coord2(:,1), embed_coord2(:,2),'.', 'color', 0.5*ones(1,3))
plot(embed_coord(SubDomain_edges,1), embed_coord(SubDomain_edges,2),'or')

%%

make_fig;
plot3(data(:,1), data(:,2), data(:,3),'.', 'color', 0.9*ones(1,3))
hold on
plot3(data2(:,1), data2(:,2), data2(:,3),'.', 'color', 0.5*ones(1,3))
% plot3(data(SubDomain_edges,1), data(SubDomain_edges,2),data(SubDomain_edges,3),'or')
xlabel('S')
ylabel('E')
zlabel('E:S')
print('data_with_subdomain', fmt, res);


%%

make_fig;
plot3(data(:,1), data(:,2), data(:,3),'.', 'color', 0.9*ones(1,3))
hold on
plot3(data2(:,1), data2(:,2), data2(:,3),'.', 'color', 0.5*ones(1,3))
plot3(data(SubDomain_edges,1), data(SubDomain_edges,2),data(SubDomain_edges,3),'or')
xlabel('S')
ylabel('E')
zlabel('E:S')
print('data_with_edge', fmt, res);

%%
idx_edge = SubDomain_edges(50);

K = 50;
% idx_nn = knnsearch(data2,data(idx_edge,:), 'K', K);
idx_nn = knnsearch(embed_coord2,embed_coord(idx_edge,:), 'K', K);

make_fig;
plot3(data(:,1), data(:,2), data(:,3),'.', 'color', 0.9*ones(1,3))
hold on
plot3(data2(:,1), data2(:,2), data2(:,3),'.', 'color', 0.5*ones(1,3))
plot3(data2(idx_nn,1),data2(idx_nn,2),data2(idx_nn,3),'.r')
plot3(data(idx_edge,1), data(idx_edge,2), data(idx_edge,3),'.b', 'markersize', 20)
xlabel('S')
ylabel('E')
zlabel('E:S')
print('data_with_nn', fmt, res);

%%

delta = data(idx_edge,:) - mean(data2(idx_nn,:));

new_point = data(idx_edge,:) + 2*delta;

make_fig;
plot3(data(:,1), data(:,2), data(:,3),'.', 'color', 0.9*ones(1,3))
hold on
plot3(data2(:,1), data2(:,2), data2(:,3),'.', 'color', 0.5*ones(1,3))
plot3(data2(idx_nn,1),data2(idx_nn,2),data2(idx_nn,3),'.r')
plot3([mean(data2(idx_nn,1)) data(idx_edge,1) new_point(1)], [mean(data2(idx_nn,2)) data(idx_edge,2) new_point(2)], [mean(data2(idx_nn,3)) data(idx_edge,3) new_point(3)],'-ob')
xlabel('S')
ylabel('E')
zlabel('E:S')
print('data_with_extension', fmt, res);