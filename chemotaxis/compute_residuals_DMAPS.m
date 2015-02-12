function res = compute_residuals_DMAPS(V, eps_med_scale)
% computes the local linear regression error for each of the DMAPS
% eigenvectors as a function of the previous eigenvectors
% V are the DMAPS eigenvectors, stored in columns and ordered
% V(:,1) is assumed to be the trivial constant eigenvector
% eps_med_scale is the scale to use in the local linear regression kernel 
% the kernel will be a Gaussian with width median(distances)/eps_med_scale
% I typically take eps_med_scale = 3
% res are the residuals of each of the fitted functions
% res(1) is to be ignored, and res(2) will always be 1
% res(i) is large/close to 1 if V(:,i) parameterizes a new direction in the
% data, and res(i) is close to 0 if V(:,i) is a harmonic of a previous
% eigenvector

n = size(V, 2);

res = zeros(n,1);
res(2) = 1;

warning('off', 'MATLAB:nearlySingularMatrix');

for i=3:n
    [~, res(i)] = local_linear_regression(V(:,i), V(:, 2:i-1), eps_med_scale);
end
