function [V, D, eps, Dis] = NIV_return_dist(data, inv_c, eps, N, alpha)

[m, n] = size(data);

Dis = zeros(m);
%h = waitbar(0, 'Please wait');
for i=1:m
    %waitbar(i/m, h);
    Dis(:,i) = sum((data - repmat(data(i,:),m,1)) * inv_c(:,:,i) .* (data - repmat(data(i,:),m,1)), 2);
end
%close(h);
Dis = (Dis + Dis')/2;

if eps == 0
    eps = median(Dis(:));
end

if alpha == 0
    [V, D] = dmaps(Dis, eps, N);
else
    [V, D] = dmaps_weight(Dis, alpha, eps, N);
end

