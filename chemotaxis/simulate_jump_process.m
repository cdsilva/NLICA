clear all
close all

%% simulations

% rate of switching velocity
lambda = 10;
% speed of particles
s = 10;

% number of particles
N = 1000;
% time step
dt = 10;
% maximum simulation time
tmax = 100;
%initial probility of a particle moving to the right (vel=+1)
pinit = 0.1; 

% simulate velocity jump process
[pos, vel, time] = simulate_jumps(N, tmax, dt, lambda, s, pinit);

%% plot results

cmap = jet;
cmap = cmap(floor(linspace(1, size(cmap,1), length(time))), :);

nbins = 25;
x_hist = linspace(-s*tmax, s*tmax, nbins);
hist_all = zeros(nbins, length(time));
hist_left = zeros(nbins, length(time));
hist_right = zeros(nbins, length(time));


% histogram of all particle positions as a function of time
for i=1:length(time)
    
    [hist_all(:,i), ~] = hist(pos(:,i),x_hist);
    
    idx = (vel(:,i) == -1);
    [hist_left(:,i), ~] = hist(pos(idx,i),x_hist);
    
    idx = (vel(:,i) == 1);
    [hist_right(:,i), ~] = hist(pos(idx,i),x_hist);
    
end

figure;
set(gcf,'DefaultAxesColorOrder',cmap)
plot(x_hist, hist_all)


figure;
set(gcf,'DefaultAxesColorOrder',cmap)
plot(x_hist, hist_left)

figure;
set(gcf,'DefaultAxesColorOrder',cmap)
plot(x_hist, hist_right)
