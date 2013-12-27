function [V, D] = NLICA(data, inv_c, neigs)

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

A = exp(-Dis/(4*ep));
A = diag(1./sum(A)) * A;

[V, D] = eigs(A, neigs);
[~, I] = sort(abs(diag(D)),'descend');
V = V(:,I);
D = D(I,I);


