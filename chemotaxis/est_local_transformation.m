function [J, detJ] = est_local_transformation(Y, X)


n = size(X, 1);

J = zeros(size(X,2), size(Y,2), n);
detJ = zeros(n, 1);

K = squareform(pdist(X));
eps = median(K(:));
W = exp(-K.^2 / eps^2);

for i=1:n
    Xx = [ones(size(X,1),1) X-repmat(X(i,:), n, 1)];
    Wx = diag(W(i,:));
    for j=1:size(Y,2)
        A = inv(Xx'*Wx*Xx)*Xx'*Wx*Y(:,j);
        J(:,j,i) = A(2:end);
    end
    detJ(i) = det(J(:,:,i));
end

