clear all
close all

%%
[X, Y] = meshgrid(1:10, 1:100);
X = (X - mean(X(:)));
Y = (Y - mean(Y(:)))/100;

data = [X(:) Y(:)];

figure;

subplot(2,3,1)
plot(data(:,1),data(:,2),'.')
axis equal
title('data')

%%

[V, D] = PCA(data, 2);

subplot(2,3,2)
scatter(data(:,1),data(:,2),50, data*V(:,1),'.')
axis equal
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','PC 1');
title('PCA')

subplot(2,3,3)
scatter(data(:,1),data(:,2),50, data*V(:,2),'.')
axis equal
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','PC 2');
title('PCA')

%%
W = squareform(pdist(data)).^2;
eps = median(W(:));
[V, D] = dmaps(W, eps, 10);

subplot(2,3,4)
scatter(data(:,1),data(:,2),50, V(:,2),'.')
axis equal
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_1');
title('DMAPS')

subplot(2,3,5)
scatter(data(:,1),data(:,2),50, V(:,3),'.')
axis equal
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_2');
title('DMAPS')

subplot(2,3,6)
scatter(data(:,1),data(:,2),50, V(:,6),'.')
axis equal
cbar = colorbar('peer',gca); 
set(get(cbar,'xlabel'),'String','\phi_5');
title('DMAPS')
