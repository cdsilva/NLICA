function V_new = decouple_dmaps(V)

V_new = V;

theta_psi = atan(0.5 * (V(:, 4).' * (V(:,2).*V(:,2) - V(:,3).*V(:,3)) ) / ( V(:,4).' * (V(:,2) .* V(:,3)))) / 2;
V_new(:, 2) = V(:, 2) * cos(theta_psi) - V(:,3) * sin(theta_psi);
V_new(:, 3) = V(:, 2) * sin(theta_psi) + V(:,3) * cos(theta_psi);

