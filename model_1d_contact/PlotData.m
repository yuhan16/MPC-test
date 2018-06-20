function PlotData(xop, uop, Jop, top, fop, zop, plot_switch)

% plot state evolution
if plot_switch(1)
    figure;
    hold on
    plot(xop(1, :));
    plot(xop(2, :));
    plot(xop(3, :));
    title('position evolution');
    legend('q_1', 'q_2', 'q_3');

    figure;
    hold on
    plot(xop(4, :));
    plot(xop(5, :));
    plot(xop(6, :));
    title('veolicity evolution');
    legend('dq_1', 'dq_2', 'dq_3');
end

% plot optimal input
if plot_switch(2)
    figure;
    hold on
    plot(uop(1, :));
    plot(uop(2, :));
    plot(uop(3, :));
    title('optimal input');
    legend('u_1', 'u_2', 'u_3');
end

% plot optimal cost
if plot_switch(3)
    figure;
    hold on
    plot(Jop);
    title('optmal cost');
end

% plot time cost
if plot_switch(4)
    figure;
    hold on
    plot(top);
    title('time cost');
end

% plot contact force and sequence
if plot_switch(5)
    figure;
    hold on
    plot(fop(1,:));
    hold on
    stairs(zop(1,:));
    title('force and corresponding bin var (right side)');
    legend('force_1','bin_1')
    
    figure;
    hold on
    plot(fop(2,:));
    hold on
    stairs(zop(2,:));
    title('force and corresponding bin var (left side)');
    legend('force_2','bin_2')
end

end