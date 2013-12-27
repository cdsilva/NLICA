function [V, D] = NLICA2(data, inv_c, N, alpha)

[m, n] = size(data);

Dis = zeros(m);

h = waitbar(0, 'Please wait');
for i=1:m
    waitbar(i/m, h);
    Dis(:,i) = sum((data - repmat(data(i,:),m,1)) * inv_c(:,:,i) .* (data - repmat(data(i,:),m,1)), 2);
end
close(h);

Dis = Dis + Dis';

ep = median(median(Dis))

[V, D] = dmaps_weight(Dis, alpha, ep, N);

% A = exp(-Dis/(4*ep));
% D = sum(A);
% 
% A = diag(D.^(-alpha)) * A * diag(D.^(-alpha));
% D = sum(A);
% 
% A = diag(1./D) * A;
% 
% [V, D] = eigs(A, N);
% [~, I] = sort(abs(diag(D)),'descend');
% V = V(:,I);
% D = D(I,I);
