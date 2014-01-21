function Y = DriftFn(t, X, k)

r = rates(X, k);

Y = zeros(5, 1);

Y(1) = -r(1) + r(2);
Y(2) = -r(1) + r(2);
Y(3) = r(1) - r(2) - r(3);
Y(4) = -r(3);
Y(5) = r(3);