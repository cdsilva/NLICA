function [V2, D2] = reorder_DMAPS_noharmonics(V, D)

D2 = D;
n = size(D, 2);

for i=3:n
    [~, res] = local_linear_regression(V(:,i), V(:, 2:i-1));
    D2(i,i) = D(i,i) * res;
end

[~, I] = sort(abs(diag(D2)), 'descend');

V2 = V(:,I);
D2 = D2(I,I);