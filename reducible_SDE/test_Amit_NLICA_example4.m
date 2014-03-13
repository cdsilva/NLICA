clear all
close all

%% define parameters

dim = 2;

% initial conditions
% data0 = [0; 0];

% tmax = 2000;

dt = 0.005;
% nsteps_per_step = 100;
% 
% s1 = 100*(dt^2);
% a = 0;
% 
% drift = @(t, x) a*ones(dim, 1);
% diffn = @(t, x) s1*eye(dim);
% 
% SDE = sde(drift, diffn, 'StartState', data0);
% nPeriods = ceil(tmax / dt);
%[data, t, Z] = SDE.simulate(nPeriods, 'DeltaTime', dt, 'NSTEPS', nsteps_per_step);

npoints = 1e6;

data = randn(npoints, dim)*dt;
data(:,1) = cumsum(data(:,1));
data(:,2) = cumsum(data(:,2));

for i=1:size(data, 1)
    for j=1:2
        if data(i, j) < 0
            data(i:end, j) = -data(i:end, j);
        elseif data(i, j) > 1
            data(i:end, j) = 2-data(i:end, j);
        end
    end
end


u = data(:,1);
v = data(:,2);

%%
f = @(x) [x(:,1)+x(:,2).^3 x(:,2)-x(:,1).^3];

fdata = f(data);

x = fdata(:,1);
y = fdata(:,2);

figure; plot(x, y, '.')

%%

npoints_per_cov = 100;
stride = 500;
[inv_c, new_data, ranks] = covariances(fdata, npoints_per_cov, 2, stride);

%%

if size(new_data, 1) > 4000
    disp('ERROR')
    return
end

neigs = 10;
%eps = 2*0.005;
eps = 0;
[V, D, ~] = NIV(new_data, inv_c, eps, neigs, 0);
V = decouple_dmaps(V);

figure;
scatter(new_data(:, 1), new_data(:,2),200,V(:,2),'.')
xlabel('x_1')
ylabel('x_2')

figure;
scatter(u(1:stride:end), v(1:stride:end),200,V(:,2),'.')
xlabel('x_1')
ylabel('x_2')

figure; 
plot(V(:,2),V(:,3),'.')
xlabel('\phi_2')
ylabel('\phi_3')

%%

ind = 17000;

figure; 
plot(x, y, '.')
hold on
plot(x(ind:ind+2*npoints_per_cov), y(ind:ind+2*npoints_per_cov), '.r')

figure; 
plot(u, v, '.')
hold on
plot(u(ind:ind+2*npoints_per_cov), v(ind:ind+2*npoints_per_cov), '.r')
