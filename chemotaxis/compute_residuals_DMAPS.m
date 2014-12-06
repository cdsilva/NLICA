function res = compute_residuals_DMAPS(V, eps_med_scale)

n = size(V, 2);

res = zeros(n,1);
res(2) = 1;

for i=3:n
    [~, res(i)] = local_linear_regression(V(:,i), V(:, 2:i-1), eps_med_scale);
end
