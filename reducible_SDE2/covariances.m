function inv_c = covariances(data, knn, D)

% Find neighborhoods
[M, N] = size(data);

m = ceil(M);

% compute covariances
inv_c = zeros(N, N, m);

for i = 1:M
    lind = max(i-knn, 1);
    rind = min(i+knn, M);
    neighbors = data(lind:rind,:);

    c = cov(neighbors);

    [U S V] = svd(c);
    %c_hat = U(:,1:D) * S(1:D,1:D) * V(:,1:D)';
    %inv_c(:,:,i) = pinv(c_hat);
    inv_c(:,:,i) = U(:,1:D) * diag(1./diag(S(1:D,1:D))) * V(:,1:D)';
end

