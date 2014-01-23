function res = dist_corr(X, Y)

n = length(X);

a = squareform(pdist(X));
b = squareform(pdist(Y));

A = a - repmat(mean(a, 1), n, 1) - repmat(mean(a, 2), 1, n) + mean(a(:));
B = b - repmat(mean(b, 1), n, 1) - repmat(mean(b, 2), 1, n) + mean(b(:));

dCov = sqrt(sum(sum(A.*B))/(n^2));
dVarX = sqrt(sum(sum(A.*A))/(n^2));
dVarY = sqrt(sum(sum(B.*B))/(n^2));

res = dCov / sqrt(dVarX*dVarY);
