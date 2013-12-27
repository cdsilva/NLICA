

%%%  plots in the manner of the LP example for one of the 30 coordinates (this is the last one, z-cord of atom #10)

load train_function_paper_test;
load  Train_EstP_paper_test;
load  deltas_paper_test;

 for i=1:11
    figure; plot(train_function_paper_test(:,30), 'b')
    hold on;
    t = Train_EstP_paper_test{i};
    plot(t(:,30), 'k');
    d = deltas_paper_test{i};
    figure; plot(d(:,30), 'r');
    axis([0 4500 -0.5 0.5]);
end;