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
    x1 = d1{1};
    x2 = d2{1};
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
    u1 = d1{2};
    u2 = d2{2};
    err = abs(u1(1, :) - u2(1, :));
    figure;
    hold on
    plot(u1(1, :));
    plot(u2(1, :));
    plot(err);
    title(['input u_1 and abs error (', m1, ' vs ', m2, ')']);
    legend(m1, m2, 'error');
    fprintf('average abs error of u_1: %f\n', mean(err));
end

if plot_switch(3)
    J1 = d1{3};
    J2 = d2{3};
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
    t1 = d1{4};
    t2 = d2{4};
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
    c1 = d1{5};
    c2 = d2{5};
    err = abs(c1(1, :) - c2(1, :));
    figure;
    hold on
    plot(c1(1, :));
    plot(c2(1, :));
    plot(err, 'x');
    title(['right contact force and abs error (', m1, ' vs ', m2, ')']);
    legend(m1, m2, 'error');
    fprintf('average abs error for cn_1: %f\n', mean(err));
end
    
end