function [t, n] = Gillespie3(n0, tmax)
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

t = 0;
n = n0;

while t < tmax
    
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
    t = t + dt;
    
    cdf = cumsum(rate)/rtotal;
    k = rand;
    ind = find(k < cdf, 1, 'first');
    
    n = n + stoich_matrix(ind,:);
end

    
