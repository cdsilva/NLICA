% this function siulates the velocity jump process from  "The diffusion
% limit of transport equations deroived from velocity-jump processes," 
% Hillen and Othmer, SIAM J. Applied Math, 2000
% N is the number of particles
% tmax is the time to simulate
% dt is the time step of simulation
% lambda is the rate of switching of left/right velocities
% s is the speed (constant)
% pos stores the positions of all the particles at each timestep (N x nsteps matrix)
% vel stores the velocities (+/-1) of all the particles at each timestep (N x nsteps matrix)
% time stores the times at each timestep
% pinit is the initial probability of a particle having velocity +1
function [pos, vel, time] = simulate_jumps_unifdist(N, tmax, dt, lambda, s, pinit, dist_min, dist_max)

%% arrays to store times, positions, and velocities
time = 0:dt:tmax;
pos = zeros(N, length(time));
vel = zeros(N, length(time));

%% keep track of current time, position, and velocity during iterations

% start time at 0
ti = 0;
% start position at uniform distribution
posi = unifrnd(dist_min, dist_max, N,1);
% intial velocities are  +/-1 with probability pinit of +1
veli =  2 * binornd(1,pinit, [N 1]) - 1;

% store initial configurations
pos(:,1) = posi;
vel(:,1) = veli;

%% simulate/update positions of particles, change velocity of random particle every time an "event" occurs

% time of next event (time between events is exponentially distributed)
tevent = ti + exprnd(1/lambda);

% update position of particles at each timestep
for i=2:length(time)
    % if event occurs before next timestep, update pos/vel and draw
    % another event
    while time(i) > tevent
        % update positions
        posi = posi + veli * s * (tevent - ti);
        % switch velocity of random particle
        switch_idx = randi(N, 1);
        veli(switch_idx) = -veli(switch_idx);
        % update current time of simulation
        ti = tevent;
        % draw new event
        tevent = ti + exprnd(1/lambda);
    end
    % update positions/velocities to current timestep
    posi = posi + veli * s * (time(i) - ti);
    pos(:,i) = posi;
    vel(:,i) = veli;
    % update current time
    ti = time(i);
end 