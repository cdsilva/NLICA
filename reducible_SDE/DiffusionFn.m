function Y = DiffusionFn(t, X, k)

r = rates(X, k);

Y = zeros(5, 3);

Y(1,:) = [-sqrt(r(1)) sqrt(r(2)) 0];
Y(2,:) = [-sqrt(r(1)) sqrt(r(2)) 0];
Y(3,:) = [sqrt(r(1)) -sqrt(r(2)) -sqrt(r(3))];
Y(4,:) = [0 0 -sqrt(r(3))];
Y(5,:) = [0 0 sqrt(r(3))];