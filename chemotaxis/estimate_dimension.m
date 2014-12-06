function d = estimate_dimension(V, D, res_tol, lambda_tol)

numD = find(abs(diag(D)) > lambda_tol, 1, 'last');

d = 1;
for i=3:numD
    [~, res] = local_linear_regression(V(:,i), V(:, 2:i-1));
    if res > res_tol^2
        d = d + 1;
    end
end


