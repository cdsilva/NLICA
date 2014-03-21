clear all
close all

%% define parameters

dim1 = 2;
dim2 = 3;

s1 = 0.001;
s2 = 0.1;

L1 = 1;
L2 = 1;

npoints = 1e5;
subsample_int = 5;

theta = randn(npoints, dim1)*s1;
theta = cumsum(theta);
for i=1:npoints
    for j=1:dim1
        if theta(i, j) < 0
            theta(i:end, j) = -theta(i:end, j);
        elseif theta(i, j) > L1
            theta(i:end, j) = 2*L1-theta(i:end, j);
        end
    end
end

xi = randn(npoints, dim2)*s2;
xi = cumsum(xi);
for i=1:npoints
    for j=1:dim2
        if xi(i, j) < 0
            xi(i:end, j) = -xi(i:end, j);
        elseif xi(i, j) > L2
            xi(i:end, j) = 2*L2-xi(i:end, j);
        end
    end
end
xi = xi - L2/2;

theta = theta(1:subsample_int:end, :);
xi = xi(1:subsample_int:end, :);


figure;
plot(theta(:,1), theta(:,2),'.')

figure;
plot3(xi(:,1),xi(:,2),xi(:,3),'.')

figure;
plot(theta(:,1)+xi(:,1),'-r')
hold on
plot(smooth(theta(:,1)+xi(:,1), 401), '-g')
plot(theta(:,1))

return
%%
alpha = 4;
sin_scale = 0.2;
f = @(x) [x(:,1)+x(:,2).^3 x(:,2)-x(:,1).^3 sin_scale*sin(alpha*(x(:,1)+x(:,2).^3))+sin_scale*sin(alpha*(x(:,2)-x(:,1).^3))];


x = f(theta) + xi;

figure; plot3(x(:,1), x(:,2), x(:,3), '.')

return
%%

npoints_per_cov = 20;
stride = 250;
[inv_c, new_data, ranks] = covariances(x, npoints_per_cov, 2, stride);

%%

if size(new_data, 1) > 4000
    disp('ERROR')
    return
end

neigs = 10;
%eps = 2*0.005;
eps = 0;
[V, D, ~] = NIV(new_data, inv_c, eps, neigs, 1/2);
V = decouple_dmaps(V);

% figure;
% scatter(new_data(:, 1), new_data(:,2),200,V(:,2),'.')
% xlabel('x_1')
% ylabel('x_2')

figure;
scatter(u(1:stride:end), v(1:stride:end),200,V(:,2),'.')
xlabel('x_1')
ylabel('x_2')

figure; 
plot(V(:,2),V(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')

%%

W = squareform(pdist(new_data)).^2;
eps2 = median(W(:));
[V2, D2] = dmaps(W, eps2, neigs);
V2 = decouple_dmaps(V2);

figure;
scatter(u(1:stride:end), v(1:stride:end),200,V2(:,2),'.')
xlabel('x_1')
ylabel('x_2')

figure; 
plot(V2(:,2),V2(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')

%%

ind = 17000;

figure; 
plot(x, y, '.')
hold on
plot3(x(ind:ind+2*npoints_per_cov), y(ind:ind+2*npoints_per_cov), z(ind:ind+2*npoints_per_cov), '.r')

figure; 
plot(u, v, '.')
hold on
plot(u(ind:ind+2*npoints_per_cov), v(ind:ind+2*npoints_per_cov), '.r')
