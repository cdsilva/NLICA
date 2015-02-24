function [fx, res] = local_linear_regression(y, X, eps_med_scale)

n = size(X, 1);

K = squareform(pdist(X));
eps = median(K(:))/eps_med_scale;
W = exp(-K.^2 / eps^2);

L = zeros(n);
for i=1:n
    Xx = [ones(size(X,1),1) X-repmat(X(i,:), n, 1)];
%     Wx = diag(W(i,:));
%     A = (Xx'*Wx*Xx)\(Xx'*Wx);
    Xx2 = Xx'.*repmat(W(i,:), size(Xx, 2), 1);
    A = (Xx2*Xx)\Xx2;
    L(i,:) = A(1,:);
end

fx = L*y;
res = sqrt(mean(((y-fx)./(1-diag(L))).^2)) / std(y);