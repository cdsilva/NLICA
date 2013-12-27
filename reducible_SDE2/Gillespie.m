function [t, n_all] = Gillespie(n0, tmax, step_per_it)
% this function performs the Gillespie algorithm 
% n0 is a row-vector of initial conditions (MUST be only integers)
% tmax is the maximum time to integration
% step_per_it is the number of Gillespie steps to take between each saved
% data point
% t is a column vector with the times sampled by the algorithm
% n_all is a matrix with each row corresponding to the concentrations at
% the time in t


if nargin < 3
    step_per_it = 1;
end

parameters;

nsize = 2^8;
n_all = zeros(nsize, length(n0));
n_all(1,:) = n0;
n = n_all(1,:);
steps = 1;
t = zeros(nsize, 1);

while t(steps) < tmax
    steps = steps + 1;
    if steps > nsize
        n_all = [n_all; zeros(size(n_all))];
        t = [t; zeros(size(t))];
        nsize = 2 * nsize;
    end
    
    t(steps) = t(steps-1);
    
    for i=1:step_per_it
        %n = n_all(steps-1,:);
    
        rate = calc_rate(n, k);
    
        rtotal = sum(rate);
    
        dt = -log(rand) / rtotal;
        t(steps) = t(steps) + dt;
    
        cdf = cumsum(rate)/rtotal;
        k = rand;
        ind = find(k < cdf, 1, 'first');
        dn = stoich_matrix(ind,:);
        n = n + dn;
    end
    
    n_all(steps,:) = n;
end

n_all = n_all(1:steps,:);
t = t(1:steps);
    
    
