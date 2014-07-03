function [inv_c, new_data, ranks] = covariances2(data, D)

% Find neighborhoods
[M, N, P] = size(data);

% compute covariances
inv_c = zeros(N, N, P);

new_data = zeros(P, N);
ranks = zeros(P, 1);

center_idx = floor(M/2);

for i = 1:P    
    c = cov(data(:,:,i));
    
    [U, S, V] = svd(c);
    %c_hat = U(:,1:D) * S(1:D,1:D) * V(:,1:D)';
    %inv_c(:,:,idx) = pinv(c_hat);
    inv_c(:,:,i) = U(:,1:D) * diag(1./diag(S(1:D,1:D))) * V(:,1:D)';
    new_data(i,:) = data(center_idx, :, i);
    ranks(i) = rank(c);
end

