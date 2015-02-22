clear all
close all

%%

scale_vec = repmat(2:10, 1 ,3);

lambda_ratio = zeros(size(scale_vec));

n = 500;

figure;
for i=1:length(scale_vec)
    
    scale = scale_vec(i);
    
    
    data = rand(n, 2);
    data(:,1) = data(:,1) * scale;
    
    %     figure;
    %     scatter(data(:,1),data(:,2))
    %     axis equal
    
    W = squareform(pdist(data)).^2;
    eps = median(W(:))/9;
    
    [V, D] = dmaps(W, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
    %     figure;
    %     make_colored_bars(diag(D(2:end, 2:end)), res(2:end))
    %     xlabel('k')
    %     ylabel('\mu_k')
    
    %     idx = find(res > 0.2);
    %     D(idx, idx);
    [~, idx] = sort(res, 'descend');
    idx = idx(1:2);
    if idx(1) > idx(2)
        tmp = idx(2);
        idx(2) = idx(1);
        idx(1) = tmp;
    end
    
    subplot(length(scale_vec), 2, 2*i-1);
    plot(data(:,1),V(:,idx(1)),'.')
    
    subplot(length(scale_vec), 2, 2*i);
    plot(data(:,2),V(:,idx(2)),'.')
    
    lambda_ratio(i) = 1/sqrt(log(D(idx(1), idx(1)))/log(D(idx(2), idx(2))));
    
end

figure;
plot(scale_vec, lambda_ratio,'o');
hold on

%%
scale_vec = repmat(2:10, 1 ,3);

lambda_ratio = zeros(size(scale_vec));

n = 500;

figure;
for i=1:length(scale_vec)
    
    scale = scale_vec(i);
    
    
    data = rand(n, 2);
    data(:,1) = data(:,1) * scale;
    
    %     figure;
    %     scatter(data(:,1),data(:,2))
    %     axis equal
    
    W = squareform(pdist(data)).^2;
    eps = median(W(:))/9;
    
    [V, D] = dmaps(W, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
    %     figure;
    %     make_colored_bars(diag(D(2:end, 2:end)), res(2:end))
    %     xlabel('k')
    %     ylabel('\mu_k')
    
    %     idx = find(res > 0.2);
    %     D(idx, idx);
    [~, idx] = sort(res, 'descend');
    idx = idx(1:2);
    if idx(1) > idx(2)
        tmp = idx(2);
        idx(2) = idx(1);
        idx(1) = tmp;
    end
    
    subplot(length(scale_vec), 2, 2*i-1);
    plot(data(:,1),V(:,idx(1)),'.')
    
    subplot(length(scale_vec), 2, 2*i);
    plot(data(:,2),V(:,idx(2)),'.')
    
    lambda_ratio(i) = 1/sqrt(log(D(idx(1), idx(1)))/log(D(idx(2), idx(2))));
%     lambda_ratio(i) = 1/sqrt(sum(log(diag(D(idx(1):idx(2)-1, idx(1):idx(2)-1))))/log(D(idx(2), idx(2))));
%     lambda_ratio(i) = D(idx(1), idx(1))/D(idx(2), idx(2))*2;
%     lambda_ratio(i) = sum(diag(D(idx(1):idx(2)-1, idx(1):idx(2)-1)))/D(idx(2), idx(2))*2;
    
    
end

%%
figure;
plot(scale_vec, lambda_ratio,'o');
hold on
plot(scale_vec, scale_vec)

%%

scale_vec = 3:6;

lambda_ratio = zeros(size(scale_vec));

figure;
for i=1:length(scale_vec)
    
    scale = scale_vec(i);
    
    data = rand(n, 2);
    
    data(:,1) = data(:,1).^scale_vec(i);
    
%     data(1:round(n/scale_vec(i)),1) = data(1:round(n/scale_vec(i)),1) * 0.5;
%     data(round(n/scale_vec(i))+1:end,1) = data(round(n/scale_vec(i))+1:end,1) * 0.5 + 0.5;
%     data(:,2) = data(:,2) / 2;
    
%         figure;
%         scatter(data(:,1),data(:,2))
%         axis equal
        
    
    W = squareform(pdist(data)).^2;
    eps = median(W(:))/9;
    
    [V, D] = dmaps(W, eps, 15);
    
    eps_med_scale = 3;
    res = compute_residuals_DMAPS(V, eps_med_scale);
    
        figure;
        make_colored_bars(diag(D(2:end, 2:end)), res(2:end))
        xlabel('k')
        ylabel('\mu_k')

% figure;
% scatter(data(:,1),data(:,2), 50, V(:,2),'.')
    
        idx = find(res > 0.2);
    %     D(idx, idx);
%     [~, idx] = sort(res, 'descend');
%     idx = idx(1:2);
%     if idx(1) > idx(2)
%         tmp = idx(2);
%         idx(2) = idx(1);
%         idx(1) = tmp;
%     end
    
%     subplot(length(scale_vec), 2, 2*i-1);
%     plot(data(:,1),V(:,idx(1)),'.')
%     
%     subplot(length(scale_vec), 2, 2*i);
%     plot(data(:,2),V(:,idx(2)),'.')
    
    lambda_ratio(i) = 1/sqrt(log(D(idx(1), idx(1)))/log(D(idx(2), idx(2))));
%     lambda_ratio(i) = 1/sqrt(sum(log(diag(D(idx(1):idx(2)-1, idx(1):idx(2)-1))))/log(D(idx(2), idx(2))));
%     lambda_ratio(i) = D(idx(1), idx(1))/D(idx(2), idx(2))*2;
%     lambda_ratio(i) = sum(diag(D(idx(1):idx(2)-1, idx(1):idx(2)-1)))/D(idx(2), idx(2))*2;
    
    
end

%%
figure;
plot(scale_vec, lambda_ratio,'o');
hold on
plot(scale_vec, scale_vec)
