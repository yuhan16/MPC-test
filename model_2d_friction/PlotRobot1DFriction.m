function PlotRobot1DFriction( x )

r = 0.6;
% define axis
axis([-6 6 -6 6]);
axis square
hold on
% barrier
line([-5,-5], [-6,6], 'Color','black', 'LineWidth', 1);
line([5,5], [-6,6], 'Color','black', 'LineWidth', 1);
line([-6,6], [5,5], 'Color','black', 'LineWidth', 1);
line([-6,6], [-5,-5], 'Color','black', 'LineWidth', 1);
line([-6,6], [-6;6], 'Color','black', 'LineWidth', 1);

% robot link
% base link
rectangle('Position',[x(1)-r x(2)-r 2*r 2*r], 'Curvature', [1 1]);
plot(x(1), x(2), 'kx');
line([x(1), x(1)+r], [x(2), x(2)], 'Color', 'red');

% display data
text(-4.5, 4.5, ['q1: (', num2str(x(1), '%.2f'), ',' num2str(x(2), '%.2f'), ')']);

end