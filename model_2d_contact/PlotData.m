function PlotData(data, plot_switch)
% this function aims to plot the simulation data
% instruction:
% plot switch for each component: [state, input, cost, time, force]
% value = 1: plot, value = 0: not plot
% Notice: 
%  -- only for state and input plot, [0, 1, ..., 6] are allowed, 
%     which represent no plot, polt x1 & dx1 & u1, ..., plot x6 & dx6 & u6
%  -- only for contact force and binary var plot, [0, 1, ..., 4] are
%     allowed, which represents no plot, plot c1 & z1, ..., plot c4 & z4

if plot_switch(1)
    seq = plot_switch(1);
    x1 = data.xop(seq, :);
    x2 = data.xop(seq+6, :);
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
    f = data.fop(seq, :);
    z = data.zop(seq, :);
    figure;
    hold on
    plot(f);
    plot(z);
    title(['contact force and bin var for joint ', num2str(seq)]);
    legend(['c_', num2str(seq)], ['z_', num2str(seq)]);
end
    
end