function PlotData(data, plot_switch)
% this function aims to plot the simulation data
% instruction:
% plot switch for each component: [state, input, cost, time, force]
% value = 1: plot, value = 0: not plot
% Notice: 
%  -- only for state and input plot, [0, 1, 2] are allowed, 
%     which represent no plot, polt x1 & dx1 & u1, plot x2 & dx2 & u2
%  -- only for contact force and binary var plot, and friction force plot, 
%     [0, 1, ..., 3] are allowed, which represents no plot, plot c1 & z1, and
%     f11, f12, ..., plot c3 & z3, and f31, f32

if plot_switch(1)
    seq = plot_switch(1);
    x1 = data.xop(seq, :);
    x2 = data.xop(seq+2, :);
    figure;
    hold on
    plot(x1);
    plot(x2);
    title(['q_', num2str(seq), 'and dq_', num2str(seq), ' evolution']);
    legend(['q_', num2str(seq)], ['dq_', num2str(seq)]);
end

if plot_switch(2)
    seq = plot_switch(1);
    u = data.uop(seq, :);
    figure;
    hold on
    plot(u);
    title(['input u_',  num2str(seq)]);
end

if plot_switch(3)
    J = data.Jop;
    figure;
    hold on
    plot(J);
    title('optimal cost');
end

if plot_switch(4)
    t = data.top;
    figure;
    hold on
    plot(t);
    title('solving time');
end

if plot_switch(5)
    seq = plot_switch(5);
    f = data.cop(seq, :);
    z = data.zop(seq, :);
    figure;
    hold on
    plot(f);
    plot(z);
    title(['contact force and bin var for joint ', num2str(seq)]);
    legend(['c_', num2str(seq)], ['z_', num2str(seq)]);
end

if plot_switch(6)
    seq = plot_switch(5);
    f1 = data.fop(2*seq-1, :);
    f2 = data.fop(2*seq, :);
    figure;
    hold on
    plot(f1);
    plot(f2);
    title(['friction force on ', num2str(seq), ' contact surface']);
    legend('f_1', 'f_2');
end
    
end