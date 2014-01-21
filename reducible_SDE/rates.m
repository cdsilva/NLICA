function r = rates(X, k)

r = zeros(3, 1);
r(1) = k(1)*X(1)* X(2);
r(2) = k(2)*X(3);
r(3) = k(3)*X(3)*X(4);
