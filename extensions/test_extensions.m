clear all
close all

%%

fmt = '-djpeg';
res = '-r600';

load data_for_elio.mat

load EDGES.mat

[m, n] = size(data);

data2 = data(subregion_idx,:);
embed_coord2 = embed_coord(subregion_idx,:);

outside_idx = setdiff((1:m)', subregion_idx);

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
% print('data_with_subdomain', fmt, res);


%%

make_fig;
plot3(data(:,1), data(:,2), data(:,3),'.', 'color', 0.9*ones(1,3))
hold on
plot3(data2(:,1), data2(:,2), data2(:,3),'.', 'color', 0.5*ones(1,3))
plot3(data(SubDomain_edges,1), data(SubDomain_edges,2),data(SubDomain_edges,3),'or')
xlabel('S')
ylabel('E')
zlabel('E:S')
% print('data_with_edge', fmt, res);

%%
idx_edge = SubDomain_edges(95);

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
% print('data_with_nn', fmt, res);

axis_lim = 1e4;

make_fig;
plot3(data(:,1), data(:,2), data(:,3),'.', 'color', 0.9*ones(1,3))
hold on
plot3(data2(:,1), data2(:,2), data2(:,3),'.', 'color', 0.5*ones(1,3))
plot3(data2(idx_nn,1),data2(idx_nn,2),data2(idx_nn,3),'.r')
plot3(data(idx_edge,1), data(idx_edge,2), data(idx_edge,3),'.b', 'markersize', 20)
xlabel('S')
ylabel('E')
zlabel('E:S')
axis([data(idx_edge,1)-axis_lim data(idx_edge,1)+axis_lim data(idx_edge,2)-axis_lim data(idx_edge,2)+axis_lim data(idx_edge,3)-axis_lim data(idx_edge,3)+axis_lim])

%%

[V, D] = PCA(data2(idx_nn,:), 2);

make_fig;
plot3(data(:,1), data(:,2), data(:,3),'.', 'color', 0.6*ones(1,3))
hold on
plot3(data2(:,1), data2(:,2), data2(:,3),'.', 'color', 0.3*ones(1,3))
hold on
plot3(data2(idx_nn,1),data2(idx_nn,2),data2(idx_nn,3),'.r')
plot3(data(idx_edge,1), data(idx_edge,2), data(idx_edge,3),'.b', 'markersize', 20)
plane_scale = 1e4;
corners = plane_scale*[V(:,1)+V(:,2) V(:,1)-V(:,2) -V(:,1)-V(:,2) -V(:,1)+V(:,2)]; 
fill3(corners(1,:)+data(idx_edge, 1), corners(2,:)+data(idx_edge, 2), corners(3,:)+data(idx_edge, 3), 'r', 'facealpha', 0.1)
xlabel('S')
ylabel('E')
zlabel('E:S')
axis([data(idx_edge,1)-axis_lim data(idx_edge,1)+axis_lim data(idx_edge,2)-axis_lim data(idx_edge,2)+axis_lim data(idx_edge,3)-axis_lim data(idx_edge,3)+axis_lim])


data_for_PCA = data - repmat(mean(data2(idx_nn,:)), m, 1);

denoised_data_nn = (data2(idx_nn,:) - repmat(mean(data2(idx_nn,:)), K, 1)) * V * V' + repmat(mean(data2(idx_nn,:)), K, 1);
denoised_data = (data - repmat(mean(data2(idx_nn,:)), m, 1)) * V * V' + repmat(mean(data2(idx_nn,:)), m, 1);


delta = data(idx_edge,:) - mean(denoised_data_nn);

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
axis([data(idx_edge,1)-axis_lim data(idx_edge,1)+axis_lim data(idx_edge,2)-axis_lim data(idx_edge,2)+axis_lim data(idx_edge,3)-axis_lim data(idx_edge,3)+axis_lim])

make_fig;
plot(data*V(:,1), data*V(:,2),'.', 'color', 0.9*ones(1,3))
hold on
plot(data2*V(:,1), data2*V(:,2),'.', 'color', 0.5*ones(1,3))
plot(data2(idx_nn,:)*V(:,1),data2(idx_nn,:)*V(:,2),'.r')
    i = idx_edge;
        delta = data(i,:) - mean(denoised_data_nn);
        new_point = data(i,:) + 2*delta;
        plot([data(i,:)*V(:,1) new_point*V(:,1)], [data(i,:)*V(:,2) new_point*V(:,2)],'-ob')
axis([data(idx_edge,:)*V(:,1)-axis_lim data(idx_edge,:)*V(:,1)+axis_lim data(idx_edge,:)*V(:,2)-axis_lim data(idx_edge,:)*V(:,2)+axis_lim ])


make_fig;
plot3(data(:,1), data(:,2), data(:,3),'.', 'color', 0.9*ones(1,3))
hold on
plot3(data2(:,1), data2(:,2), data2(:,3),'.', 'color', 0.5*ones(1,3))
plot3(data2(idx_nn,1),data2(idx_nn,2),data2(idx_nn,3),'.r')
for j = 1:length(SubDomain_edges)
    i = SubDomain_edges(j);
    if norm(data(i, :) - data(idx_edge, :)) < 4000
        delta = data(i,:) - mean(denoised_data_nn);
        new_point = data(i,:) + 2*delta;
        plot3([data(i,1) new_point(1)], [data(i,2) new_point(2)], [data(i,3) new_point(3)],'-ob')
    end
end
axis([data(idx_edge,1)-axis_lim data(idx_edge,1)+axis_lim data(idx_edge,2)-axis_lim data(idx_edge,2)+axis_lim data(idx_edge,3)-axis_lim data(idx_edge,3)+axis_lim])

make_fig;
plot(data*V(:,1), data*V(:,2),'.', 'color', 0.9*ones(1,3))
hold on
plot(data2*V(:,1), data2*V(:,2),'.', 'color', 0.5*ones(1,3))
plot(data2(idx_nn,:)*V(:,1),data2(idx_nn,:)*V(:,2),'.r')
for j = 1:length(SubDomain_edges)
    i = SubDomain_edges(j);
    if norm(data(i, :) - data(idx_edge, :)) < 4000
        delta = data(i,:) - mean(denoised_data_nn);
        new_point = data(i,:) + 2*delta;
        plot([data(i,:)*V(:,1) new_point*V(:,1)], [data(i,:)*V(:,2) new_point*V(:,2)],'-ob')
    end
end
axis([data(idx_edge,:)*V(:,1)-axis_lim data(idx_edge,:)*V(:,1)+axis_lim data(idx_edge,:)*V(:,2)-axis_lim data(idx_edge,:)*V(:,2)+axis_lim ])


err = sum((data - denoised_data).^2, 2);
dist_from_edge = pdist2(data, data(idx_edge, :));
err2 = sum((data2(idx_nn,:) - denoised_data_nn).^2, 2);
dist_from_edge2 = pdist2(data2(idx_nn,:), data(idx_edge, :));
figure;
semilogy(dist_from_edge(outside_idx), err(outside_idx), '.')
hold on 
semilogy(dist_from_edge2, err2, '.r')

% figure;
% scatter3(data(:,1),data(:,2),data(:,3), 50, log(dist_from_edge+0.01), '.')

%%
make_fig;
plot3(data(:,1),data(:,2),data(:,3), '.', 'color', 0.9*ones(1, 3))
hold on
scatter3(data(outside_idx,1),data(outside_idx,2),data(outside_idx,3), 50, log(err(outside_idx)), '.')
plot3(data(idx_edge,1),data(idx_edge,2),data(idx_edge,3), '.k', 'markersize', 30)

make_fig;
plot(data*V(:,1),data*V(:,2), '.', 'color', 0.9*ones(1, 3))
hold on
scatter(data(outside_idx,:)*V(:,1),data(outside_idx,:)*V(:,2), 50, log(err(outside_idx)), '.')
plot(data(idx_edge,:)*V(:,1),data(idx_edge,:)*V(:,2), '.k', 'markersize', 30)

return



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
% print('data_with_extension', fmt, res);


%%
delta = data(idx_edge,:) - mean(data2(idx_nn,:));

new_point = data(idx_edge,:) + 2*delta;

make_fig;
plot3(data(:,1), data(:,2), data(:,3),'.', 'color', 0.9*ones(1,3))
hold on
plot3(data2(:,1), data2(:,2), data2(:,3),'.', 'color', 0.5*ones(1,3))
plot3(data2(idx_nn,1),data2(idx_nn,2),data2(idx_nn,3),'.r')
plot3([mean(data2(idx_nn,1)) data(idx_edge,1) new_point(1)], [mean(data2(idx_nn,2)) data(idx_edge,2) new_point(2)], [mean(data2(idx_nn,3)) data(idx_edge,3) new_point(3)],'-ob')
axis_lim = 1e4;
axis([data(idx_edge,1)-axis_lim data(idx_edge,1)+axis_lim data(idx_edge,2)-axis_lim data(idx_edge,2)+axis_lim data(idx_edge,3)-axis_lim data(idx_edge,3)+axis_lim])
xlabel('S')
ylabel('E')
zlabel('E:S')
grid on
