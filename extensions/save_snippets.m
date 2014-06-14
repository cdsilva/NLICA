% this code generates k initial conditions, and integrates them forward and
% saves the last ncov_points points in the matrix n (and the corresponding
% times in t) so that we can estimate the local covariances
% ncov_points sets how many points in the trajectory we want to save for
% each run
% k sets the number of runs (points on the manifold) we want
% use_ODE sets whether we want to integrate the ODE to find points on the
% slow manifold (faster) or use Gillespie to find points on the slow
% manifold

%% parallel operations
matlabpool open 12

%% load parameters and define constants

% the length of the short bursts to estimate the covariances (in time units)
snippet_len = 0.2;

% number of bursts 
cloud_size = 20;

% number of points to take on the manifold
k = 3000;

% rate constants, etc., defined in this file
parameters;

% time to integrate to get points on manifold
tmax = 10;

% tolerance for ode integrator
ode_opts.tol = 1e-9;

%% define set of initial conditions

% generate random sample points in 6-dimensional space
n0_all = rand(50*k, 6);
n0_all(:,1) = n0_all(:,1) * ST;
n0_all(:,2) = n0_all(:,2) * ET;
n0_all(:,3) = n0_all(:,3) * min(ST, ET);
n0_all(:,4) = n0_all(:,4) * min(ST, ET);
n0_all(:,5) = n0_all(:,5) * DT;
n0_all(:,6) = n0_all(:,6) * FT;

% find points that satisfy constraints of problem (all concentrations must
% be >=0)
ind = find(n0_all(:,1)+n0_all(:,3)+n0_all(:,4)+n0_all(:,5) < ST &...
    n0_all(:,2)+n0_all(:,3)+n0_all(:,4)+n0_all(:,6) < ET & ...
    n0_all(:,5) < DT & ...
    n0_all(:,6) < FT, k, 'first');

k = length(ind);
n0_all = n0_all(ind,:);

%% normalize    
%for j=1:6
%    n0_all(:,j) = n0_all(:,j)/max(n(:,j));
%end

%% integrate forward using Gillespie

% n will store the end points of the bursts from each initial point on the manifold
data = zeros(k, 6);
n = zeros(cloud_size, 6, k);

% integrate forward using Gillespie; store the end points of the bursts
n0_all = floor(n0_all);
parfor i=1:k
    disp(num2str(i));
    [TOUT, nOUT] = Gillespie3(n0_all(i,:),tmax)
    data(i,:) = nOUT;
    for j=1:cloud_size
        interval = 1+(j-1)*snippet_len:j*snippet_len;
        %[t(interval,i), n(interval,:,i)] = Gillespie2(floor(nOUT(end,:)),snippet_len);
        [tmp_t, tmp_n] = Gillespie3(nOUT(end,:),snippet_len);
        n(j,:,i) = tmp_n;
    end
end

save('snippet_data.mat');

matlabpool close
