function [t, n_all] = Gillespie2(n0, nsteps)
% this function performs the Gillespie algorithm for the set of reactions
% defined in "Stability and Stabilization of Constrained Runs"
% n0 is a row-vector of initial conditions (MUST be only integers)
% nsteps is the number of steps to take
% t is a column vector with the times sampled by the algorithm
% n_all is a matrix with each row corresponding to the concentrations at
% the time in t

parameters;

stoich_matrix = [-1 -1  1  0  0  0;
    1  1 -1  0  0  0;
    0  1 -1  0  0  0;
    -1 -1  0  1  0  0;
    1  1  0 -1  0  0;
    1  0  0 -1  0  0;
    0  0  0  0  1  0;
    0  0  0  0 -1  0;
    1  0  0  0 -1  0;
    0  0  0  0  0  1;
    0  0  0  0  0 -1;
    0  1  0  0  0 -1];

n_all = zeros(nsteps, length(n0));
t = zeros(nsteps, 1);

n_all(1,:) = n0;

for i=2:nsteps
    
    n = n_all(i-1,:);
    rate = [e1*n(1)*n(2);
        e_1*n(3);
        e2*n(3);
        b1*n(2)*n(1);
        b_1*n(4);
        b2*n(4);
        d1*(DT-n(5))*(ST-n(1)-n(3)-n(4)-n(5));
        d_1*n(5);
        d2*n(5);
        f1*(FT-n(6))*(ET-n(2)-n(3)-n(4)-n(6));
        f_1*n(6);
        f2*n(6)];
    
    rtotal = sum(rate);
    
    dt = exprnd(1/rtotal);
    t(i) = t(i-1) + dt;
    
    cdf = cumsum(rate)/rtotal;
    k = rand;
    ind = find(k < cdf, 1, 'first');
    dn = stoich_matrix(ind,:);
    
    n_all(i,:) = n_all(i-1,:) + dn;
end

    
    
