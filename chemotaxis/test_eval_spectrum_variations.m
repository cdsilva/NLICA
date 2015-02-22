function test_eval_spectrum_variations

clear all
close all

%% strips of different ratios with uniform points

ymax = 1;

for xmax = [2 4 8]
    
    spacing = 0.025*xmax;
    [X, Y] = meshgrid(0:spacing:xmax, 0:spacing:ymax);
    data = [X(:) Y(:)];
    
    
    W = squareform(pdist(data)).^2;
    eps = median(W(:));
    
    [V, D] = dmaps(W, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
    make_plot(data, V, D, res)

end

%% strips of different ratios with uniformly random points

n = 500;
ymax = 1;

rng(12321);

for xmax = [2 4 8]
    
    data = rand(n, 2);
    data(:,1) = data(:,1) * xmax;
    

    
    W = squareform(pdist(data)).^2;
    eps = median(W(:));
    
    [V, D] = dmaps(W, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);

        make_plot(data, V, D, res)

end

%% strips of different ratios with points dense in x

ymax = 1;
spacingx = 0.05;
spacingy = 0.5;

for xmax = [2 4 8]
    
    [X, Y] = meshgrid(0:spacingx:xmax, 0:spacingy:ymax);
    data = [X(:) Y(:)];

    
    W = squareform(pdist(data)).^2;
    eps = median(W(:));
    
    [V, D] = dmaps(W, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
        make_plot(data, V, D, res)

    
end

%% strips of different ratios with points dense in x; alpha = 1

alpha = 1;

ymax = 1;
spacingx = 0.05;
spacingy = 0.5;

for xmax = [2 4 8]
    
    
    [X, Y] = meshgrid(0:spacingx:xmax, 0:spacingy:ymax);
    data = [X(:) Y(:)];
    
    
    W = squareform(pdist(data)).^2;
    eps = median(W(:));
    
    [V, D] = dmaps_weight(W, alpha, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
        make_plot(data, V, D, res)

    
end

%% strips of different ratios with points dense in y

ymax = 1;
spacingx = 0.5;
spacingy = 0.05;

for xmax = [2 4 8]
    
    [X, Y] = meshgrid(0:spacingx:xmax, 0:spacingy:ymax);
    data = [X(:) Y(:)];
    

    
    W = squareform(pdist(data)).^2;
    eps = median(W(:));
    
    [V, D] = dmaps(W, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
        make_plot(data, V, D, res)

    
end

%% strips of different ratios with points dense in y; alpha = 1

alpha = 1;

ymax = 1;
spacingx = 0.5;
spacingy = 0.05;

for xmax = [2 4 8]
    
    [X, Y] = meshgrid(0:spacingx:xmax, 0:spacingy:ymax);
    data = [X(:) Y(:)];
    
   
    W = squareform(pdist(data)).^2;
    eps = median(W(:));
    
    [V, D] = dmaps_weight(W, alpha, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
        make_plot(data, V, D, res)

    
end

%% nonuniform distribution

npoints = 500;
ymax = 1;

rng(123);

ratio = 8;

for xmax = [2 4 8]
     box_scale = 0.5;
     
    nboxes = xmax/box_scale + ymax/box_scale;
    npoints_per_box = round(npoints/((ratio+1)*nboxes));
    
   
    
    data = [];
    for i=box_scale:box_scale:xmax
        for j=box_scale:box_scale:ymax
            if mod(i+j, 1) == 0
                data_tmp = rand(npoints_per_box, 2);
            else
                data_tmp = rand(ratio*npoints_per_box, 2);
            end
            data_tmp = data_tmp * box_scale;
            data_tmp(:,1) = data_tmp(:,1) + (i-box_scale);
            data_tmp(:,2) = data_tmp(:,2) + (j-box_scale);
            data = [data; data_tmp];
        end
    end
   

    W = squareform(pdist(data)).^2;
    eps = median(W(:))*2;
    
    [V, D] = dmaps(W, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
        make_plot(data, V, D, res)

    
end

%% nonuniform sampling; alpha = 1

alpha = 1;

npoints = 500;
ymax = 1;

rng(123);

ratio = 8;

for xmax = [2 4 8]
     box_scale = 0.5;
     
    nboxes = xmax/box_scale + ymax/box_scale;
    npoints_per_box = round(npoints/((ratio+1)*nboxes));
    
   
    
    data = [];
    for i=box_scale:box_scale:xmax
        for j=box_scale:box_scale:ymax
            if mod(i+j, 1) == 0
                data_tmp = rand(npoints_per_box, 2);
            else
                data_tmp = rand(ratio*npoints_per_box, 2);
            end
            data_tmp = data_tmp * box_scale;
            data_tmp(:,1) = data_tmp(:,1) + (i-box_scale);
            data_tmp(:,2) = data_tmp(:,2) + (j-box_scale);
            data = [data; data_tmp];
        end
    end
    
    W = squareform(pdist(data)).^2;
    eps = median(W(:))*2;
    
    [V, D] = dmaps_weight(W, alpha, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
    make_plot(data, V, D, res)
    
end

function make_plot(data, V, D, res)

figure;

% subplot(2,2,1)
subplot(1,2,1)
plot(data(:,1),data(:,2),'.')
axis equal

% subplot(2,2,2)
subplot(1,2,2)
make_colored_bars(diag(D(2:end, 2:end)), res(2:end))
xlabel('k')
ylabel('\mu_k')

idx = find(res > 0.4);
title(sprintf('evals: %2.2f, %2.2f', D(idx(1), idx(1)), D(idx(2), idx(2))))

% subplot(2,2,3)
% scatter(data(:,1),data(:,2),50,V(:,idx(1)),'.')
% axis equal
% 
% subplot(2,2,4)
% scatter(data(:,1),data(:,2),50,V(:,idx(2)),'.')
% axis equal
