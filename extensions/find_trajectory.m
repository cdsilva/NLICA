clear all
close all

load ../rxn_example/snippet_data.mat
load data_for_elio.mat

%% integrate forward using Gillespie

% integrate forward using Gillespie; store the end points of the bursts
n0_all = floor(n0_all);

i = 1016;
tmax_new = 1.5*tmax;

disp(num2str(i));
[TOUT, nOUT] = Gillespie(n0_all(i,:),tmax_new);

boundary_idx = find(TOUT > 10, 1,'first');

%%

figure; 
plot3(data(:,1),data(:,2),data(:,3), '.','color',0.9*ones(1,3))
hold on
plot3(data(subregion_idx,1),data(subregion_idx,2),data(subregion_idx,3), '.r')
plot3(nOUT(:,1),nOUT(:,2),nOUT(:,3), 'linewidth', 2)

%%
stride = 10;
trajectory =  nOUT(1:stride:end,:);
time = TOUT(1:stride:end);
boundary_idx = find(time > tmax, 1,'first');

save('trajectory_for_elio.mat', 'trajectory','time','boundary_idx')

%%
figure; 
plot3(data(:,1),data(:,2),data(:,3), '.','color',0.9*ones(1,3))
hold on
plot3(data(subregion_idx,1),data(subregion_idx,2),data(subregion_idx,3), '.r')
plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3), 'linewidth', 2)