

% Find neighborhoods
M = size(data, 1);
N = size(data, 2);
%knn = 5;
knn = 10;

% compute covariances
inv_c_noisy = zeros(N, N, M);
inv_c = zeros(N, N, M);
D=5;

for i = 1:M
    lind = max(i-knn, 1+100*floor((i-1)/100));
    rind = min(i+knn, 100*ceil(i/100));
    neighbors = data(lind:rind,:);

    c = cov(neighbors);
    inv_c_noisy(:,:,i) = pinv(c);

    [U S V] = svd(c);
    c_hat = U(:,1:D) * S(1:D,1:D) * V(:,1:D)';
    inv_c(:,:,i) = pinv(c_hat);
end

