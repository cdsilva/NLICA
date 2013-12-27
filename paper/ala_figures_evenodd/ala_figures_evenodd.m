clear all 
close all

fmt = '-djpeg';
res = '-r600';
set(0,'DefaultAxesFontSize',6)
set(0,'DefaultFigurePaperUnits','inches')

%% Load Data
%load '~/resesarch/Alanine/Short_Trajs_Cartesian_Coord.dat'
%load '~/resesarch/Alanine/Short_Trajs_dihedrals.dat'
load '../../Ronen/Short_Trajs_Cartesian_Coord.dat'
load '../../Ronen/Short_Trajs_dihedrals.dat'

time_indices = 10001:20000;

coor_indices2 = [14,15,16,20,21,22,26,27,28]; % only odd
coor_indices3 = [5,6,7,11,12,13,17,18,19,23,24,25,29,30,31]; % only even

%% Odd Coordinates
data = Short_Trajs_Cartesian_Coord(time_indices, coor_indices2);
    
% Compute local covariances
covariances;
    
% NLICA
NLICA;
for i=1:size(psi_mat,2)
    psi_mat(:,i) = sign(mean(psi_mat(:,i))) * psi_mat(:,i) / norm(psi_mat(:,i)) * size(psi_mat,1);
end
psi_mat_nlica_odd = psi_mat;
Psi1_nlica_odd = Psi_1;
Psi2_nlica_odd = Psi_2;

% DM
DM;   
for i=1:size(psi_mat,2)
    psi_mat(:,i) = sign(mean(psi_mat(:,i))) * psi_mat(:,i) / norm(psi_mat(:,i)) * size(psi_mat,1);
end
psi_mat_dm_odd = psi_mat;
Psi1_dm_odd = Psi_1;
Psi2_dm_odd = Psi_2;

%% Even Coordinates  
data = Short_Trajs_Cartesian_Coord(time_indices, coor_indices3);
    
% Compute local covariances
covariances;
    
% NLICA
NLICA;
for i=1:size(psi_mat,2)
    psi_mat(:,i) = sign(mean(psi_mat(:,i))) * psi_mat(:,i) / norm(psi_mat(:,i)) * size(psi_mat,1);
end
psi_mat_nlica_even = psi_mat;
Psi1_nlica_even = Psi_1;
Psi2_nlica_even = Psi_2;

% DM
DM;   
for i=1:size(psi_mat,2)
    psi_mat(:,i) = sign(mean(psi_mat(:,i))) * psi_mat(:,i) / norm(psi_mat(:,i)) * size(psi_mat,1);
end
psi_mat_dm_even = psi_mat;
Psi1_dm_even = Psi_1;
Psi2_dm_even = Psi_2;

%%
disp('DM correlation:')
corr(psi_mat_dm_even(:,1:3),psi_mat_dm_odd(:,1:3))

disp('NLICA correlation:')
corr(psi_mat_nlica_even(:,1:3),psi_mat_nlica_odd(:,1:3))

for i=1:3
    if corr(psi_mat_dm_even(:,i),psi_mat_dm_odd(:,i)) < 0
        psi_mat_dm_even(:,i) = -psi_mat_dm_even(:,i);
    end
    if corr(psi_mat_nlica_even(:,i),psi_mat_nlica_odd(:,i)) < 0
        psi_mat_nlica_even(:,i) = -psi_mat_nlica_even(:,i);
    end
end
        

figure;
set(gcf,'paperposition',[0 0 3 2])
plot(psi_mat_dm_even(:,2),psi_mat_dm_odd(:,2),'.','markersize',1)
xlabel('DM_2^{even}')
ylabel('DM_2^{odd}')
saveas(gcf,'DM_corr','fig')
print('DM_corr',fmt,res)

figure;
set(gcf,'paperposition',[0 0 3 2])
plot(psi_mat_nlica_even(:,2),psi_mat_nlica_odd(:,2),'.','markersize',1)
xlabel('\psi_2^{even}')
ylabel('\psi_2^{odd}')
saveas(gcf,'NLICA_corr','fig')
print('NLICA_corr',fmt,res)