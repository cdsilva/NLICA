function [t, n] = Gillespie3(n0, tmax)
% this function performs the Gillespie algorithm for the set of reactions
% defined in "Stability and Stabilization of Constrained Runs"
% n0 is a row-vector of initial conditions (MUST be only integers)
% tmax is the maximum time to integration
% t is a column vector with the times sampled by the algorithm
% n_all is a matrix with each row corresponding to the concentrations at
% the time in t

parameters;

t = 0;
n = n0;

while t < tmax
    
    rate = calc_rate(n, k);
    
    rtotal = sum(rate);
    
    dt = -log(rand) / rtotal;
    t = t + dt;
    
    cdf = cumsum(rate)/rtotal;
    k = rand;
    ind = find(k < cdf, 1, 'first');
    
    n = n + stoich_matrix(ind,:);
end

    
