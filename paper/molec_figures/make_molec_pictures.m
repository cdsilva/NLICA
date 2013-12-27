clear all
close all

%fmt = '-djpeg';
fmt = '-depsc';
res = '';
%res = '-r600';
set(0,'DefaultAxesFontSize',6)
set(0,'DefaultFigurePaperUnits','inches')

load ../ala_figures_allodd/LapPyr_recon.mat;

ala2 = pdbread('ala2_noH.pdb');

test_ind = 900;

true_data = reshape(true_pos(test_ind,:),3,[])';
recon_data = reshape(recon_pos(test_ind,:),3,[])';

for i=1:10
    ala2.Model.Atom(1,i).X = true_data(i,1);
    ala2.Model.Atom(1,i).Y = true_data(i,2);
    ala2.Model.Atom(1,i).Z = true_data(i,3);
    ala2.Model.Atom(1,i+10).X = recon_data(i,1);
    ala2.Model.Atom(1,i+10).Y = recon_data(i,2);
    ala2.Model.Atom(1,i+10).Z = recon_data(i,3);
end

h = molviewer(ala2);
evalrasmolscript(h,'rotate Y 135');
evalrasmolscript(h,'background white');

evalrasmolscript(h,'select *U; color black');
evalrasmolscript(h,'select *X; color red');
evalrasmolscript(h,'select all');

%saveas(gcf,sprintf('molecule%d_balls',test_ind),'fig')
print(sprintf('molecule%d_balls',test_ind),fmt,res)

%set(gcf,'paperposition',[0 0 2 2])
evalrasmolscript(h,'wireframe 10');

%saveas(gcf,sprintf('molecule%d_wire',test_ind),'fig')
print(sprintf('molecule%d_wire',test_ind),fmt,res)