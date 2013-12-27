function inv_c = calc_covar(data_bursts, D)

% Use end of snippets
[~, N, M] = size(data_bursts);

% compute covariances
%inv_c_noisy = zeros(N, N, M);
inv_c = zeros(N, N, M);

%spec = zeros(N,M);
for i = 1:M

    c = cov(data_bursts(:,:,i));
    %inv_c_noisy(:,:,i) = pinv(c);

    [U S V] = svd(c);
    c_hat = U(:,1:D) * S(1:D,1:D) * V(:,1:D)';
    inv_c(:,:,i) = pinv(c_hat);
    
    %spec(:,i) = diag(S);
end

