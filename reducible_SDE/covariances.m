function [inv_c, new_data, ranks] = covariances(data, knn, D, stride)

% Find neighborhoods
[M, N] = size(data);

m = floor(M / stride);

% compute covariances
inv_c = zeros(N, N, m);

new_data = zeros(m, N);
ranks = zeros(m, 1);

for i = 1:stride:M
    idx = ceil(i / stride);
    
    lind = max(i-knn, 1);
    rind = min(i+knn, M);
    
    neighbors = data(lind:rind,:);
    
    c = cov(neighbors);
    
    [U, S, V] = svd(c);
    %c_hat = U(:,1:D) * S(1:D,1:D) * V(:,1:D)';
    %inv_c(:,:,idx) = pinv(c_hat);
        inv_c(:,:,idx) = U(:,1:D) * diag(1./diag(S(1:D,1:D))) * V(:,1:D)';
    new_data(idx,:) = data(i,:);
    ranks(idx) = rank(c);
    
    
    %     if any(isnan(inv_c(:)))
    %         inv_c(:,:,i) = U(:,1:D-1) * diag(1./diag(S(1:D-1,1:D-1))) * V(:,1:D-1)';
    %         fprintf('reduced dim by 1 at %d', i);
    %     end
end

