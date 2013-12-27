fmt = '-depsc';
res = '-r600';
set(0,'DefaultAxesFontSize',6)
set(0,'DefaultFigurePaperUnits','inches')

%%%  plots in the manner of the LP example for one of the 30 coordinates (this is the last one, z-cord of atom #10)

load train_function_paper_test;
load  Train_EstP_paper_test;
load  deltas_paper_test;
load train_set_paper_test;

%for i=1:11
for i=8
    figure;
    set(gcf,'paperposition',[0 0 3 2.5])
    plot(train_set_paper_test(:,1),train_function_paper_test(:,30), '.b')
    hold on;
    t = Train_EstP_paper_test{i};
    plot(train_set_paper_test(:,1),t(:,30), '.k');
    xlabel('\psi_1')
    ylabel('f(x)')
    title(sprintf('Scale %d',i))
    saveas(gcf,sprintf('scale_recon%d',i),'fig')
    print(sprintf('scale_recon%d',i), fmt,res)
    
    
    d = deltas_paper_test{i};
    figure;
    set(gcf,'paperposition',[0 0 3 2.5])
    plot(train_set_paper_test(:,1),d(:,30), '.r');
    axis([-300 100 -0.5 0.5]);
    xlabel('\psi_1')
    ylabel('f^{LP}(x)-f(x)')
    title(sprintf('Scale %d',i))
    saveas(gcf,sprintf('scale_error%d',i),'fig')
    print(sprintf('scale_error%d',i), fmt,res)
end;