function [t, n_all] = Gillespie(n0, tmax)
% this function performs the Gillespie algorithm for the set of reactions
% defined in "Stability and Stabilization of Constrained Runs"
% n0 is a row-vector of initial conditions (MUST be only integers)
% tmax is the maximum time to integration
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

nsize = 2^8;
n_all = zeros(nsize, length(n0));
n_all(1,:) = n0;
steps = 1;
t = zeros(nsize, 1);

while t(steps) < tmax
    steps = steps + 1;
    if steps > nsize
        n_all = [n_all; zeros(size(n_all))];
        t = [t; zeros(size(t))];
        nsize = 2 * nsize;
    end
    
    n = n_all(steps-1,:);
    
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
    
    dt = -log(rand) / rtotal;
    t(steps) = t(steps-1) + dt;
    
    cdf = cumsum(rate)/rtotal;
    k = rand;
    ind = find(k < cdf, 1, 'first');
    dn = stoich_matrix(ind,:);
    
    n_all(steps,:) = n_all(steps-1,:) + dn;
end

n_all = n_all(1:steps,:);
t = t(1:steps);
    
    
