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
function [pos, vel, time] = simulate_jumps2(N, tmax, dt, lambda, s, pinit)

%% arrays to store times, positions, and velocities
time = 0:dt:tmax;
pos = zeros(N, length(time));
vel = zeros(N, length(time));

%% simulate/update positions of particles, change velocity of random particle every time an "event" occurs

% time of next event (time between events is exponentially distributed)
tevent = exprnd(1/lambda, 1, round(2*tmax*lambda));
while sum(tevent) < tmax
    tevent = exprnd(1/lambda, 1, 2*length(tevent));
end

idx = find(cumsum(tevent) > tmax, 1, 'first');
tevent = tevent(1:idx);

time_vec = [0 cumsum(tevent)];
vel_switch = zeros(N, length(tevent));
idx = sub2ind(size(vel_switch), randi(N, length(tevent), 1), (1:length(tevent))');
vel_switch(idx) = 1;
vel_switch = [binornd(1,pinit, [N 1]) vel_switch];
vel_vec = s * ((mod(cumsum(vel_switch, 2), 2) * 2) - 1);
pos_vec = cumsum(vel_vec(:, 1:end-1) .* repmat(time_vec(2:end)-time_vec(1:end-1), N, 1), 2);

for i=1:length(time)
    idx = find(time_vec > time(i), 1, 'first');
    pos(:,i) = pos_vec(:, idx-1);
    vel(:,i) = vel_vec(:, idx-1);
end

