function PlotError(data, methodName, plot_switch)
% this function aims to plot and compute errors between two group of data
% instruction:
% plot switch for each component: [state, input, cost, time, force]
% value = 1: plot, value = 0: not plot
% Notice: only for the state plot, [0, 1, 2, 3] are allowed, which represent 
%no plot, polt x1 & dx1, plot x2 & d2, plot x3 & dx3

if length(methodName) ~= 2
    fprintf('Wrong input argument.\n');
    return
end

m1 = methodName{1};
m2 = methodName{2};
d1 = data{1};
d2 = data{2};

if plot_switch(1)
    seq = plot_switch(1);
    x1 = d1.xop;
    x2 = d2.xop;
    err = abs(x1(seq, :) - x2(seq, :));
    figure;
    hold on
    plot(x1(seq, :));
    plot(x2(seq, :));
    plot(err);
    title(['q_', num2str(seq), ' evolution and abs error (', m1, ' vs ', m2, ')']);
    legend(m1, m2, 'error');
    fprintf('average abs error of q_%d: %f\n', seq, mean(err));
    
    err = abs(x1(seq+3, :) - x2(seq+3, :));
    figure;
    hold on
    plot(x1(seq+3, :));
    plot(x2(seq+3, :));
    plot(err);
    title(['dq_', num2str(seq), ' evolution and abs error (', m1, ' vs ', m2, ')']);
    legend(m1, m2, 'error');
    fprintf('average abs error of dq_%d: %f\n', seq, mean(err));
end

if plot_switch(2)
    seq = plot_switch(1);
    u1 = d1.uop;
    u2 = d2.uop;
    err = abs(u1(seq, :) - u2(seq, :));
    figure;
    hold on
    plot(u1(1, :));
    plot(u2(1, :));
    plot(err);
    title(['input u_', num2str(seq), ' and abs error (', m1, ' vs ', m2, ')']);
    legend(m1, m2, 'error');
    fprintf('average abs error of u_%d: %f\n', seq, mean(err));
end

if plot_switch(3)
    J1 = d1.Jop;
    J2 = d2.Jop;
    err = abs(J1 - J2);
    figure;
    hold on
    plot(J1);
    plot(J2);
    plot(err);
    title(['optimal cost and abs error (', m1, ' vs ', m2, ')']);
    legend(m1, m2, 'error');
    fprintf('average abs error of J: %f\n', mean(err));
end

if plot_switch(4)
    t1 = d1.top;
    t2 = d2.top;
    err = abs(t1 - t2);
    figure;
    hold on
    plot(t1);
    plot(t2);
    plot(err);
    title(['time and abs error (', m1, ' vs ', m2, ')']);
    legend(m1, m2, 'error');
    fprintf('average abs error for t: %f\n', mean(err));
end

if plot_switch(5)
    seq = plot_switch(5);
    c1 = d1.cop;
    c2 = d2.cop;
    err = abs(c1(seq, :) - c2(seq, :));
    figure;
    hold on
    plot(c1(seq, :));
    plot(c2(seq, :));
    plot(err, 'x');
    title(['contact force and abs error (', m1, ' vs ', m2, ')']);
    legend(m1, m2, 'error');
    fprintf('average abs error for cn_%d: %f\n', seq, mean(err));
end

if plot_switch(6)
    seq = plot_switch(6);
    f1 = d1.fop;
    f2 = d2.fop;
    err = abs(f1(2*seq-1, :) - f2(2*seq-1, :));
    figure;
    hold on
    plot(f1(2*seq, :));
    plot(f2(2*seq, :));
    plot(err, 'x');
    title(['friction force and abs error (', m1, ' vs ', m2, ')']);
    legend(m1, m2, 'error');
    fprintf('average abs error for f_%d: %f\n', seq, mean(err));
end
    
end