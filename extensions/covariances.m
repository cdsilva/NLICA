% Use end of snippets
[tmp, N, M] = size(n);

% compute covariances
inv_c_noisy = zeros(N, N, M);
inv_c_denoised = zeros(N, N, M);
D=2;

spec = zeros(N,M);
for i = 1:M

    c = cov(n(:,:,i));
    inv_c_noisy(:,:,i) = pinv(c);

    [U S V] = svd(c);
    c_hat = U(:,1:D) * S(1:D,1:D) * V(:,1:D)';
    inv_c_denoised(:,:,i) = pinv(c_hat);
    
    spec(:,i) = diag(S);
end

