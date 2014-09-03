function [c, inv_c, ranks] = covariances(data, D, deltat)

% Find neighborhoods
[M, N, P] = size(data);

% compute covariances
c = zeros(N, N, P);
inv_c = zeros(N, N, P);

ranks = zeros(P, 1);

for i = 1:P    
    c_tmp = cov(data(:,:,i))/deltat;
    
    [U, S, V] = svd(c_tmp);

    inv_c(:,:,i) = U(:,1:D) * diag(1./diag(S(1:D,1:D))) * V(:,1:D)';
    ranks(i) = rank(c_tmp);
    c(:,:,i) = c_tmp;
end

