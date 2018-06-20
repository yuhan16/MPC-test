function PlotRobot2D( x )

% define axis
axis([-6 6 -6 6]);
hold on
% barrier
line([-5,-5], [-6,6], 'Color','black', 'LineWidth', 1);
line([5,5], [-6,6], 'Color','black', 'LineWidth', 1);
line([-6,6], [5,5], 'Color','black', 'LineWidth', 1);
line([-6,6], [-5,-5], 'Color','black', 'LineWidth', 1);

% robot link
% base link
rectangle('Position',[x(1)-1 x(2)-1 2 2]);
plot(x(1), x(2), 'kx');
% end effector: q2, q3, q4; q5
a = 0.15;
rectangle('Position', [x(1)+1+x(3)-a x(2)-a 2*a 2*a], 'FaceColor', 'k');
rectangle('Position', [x(1)-a x(2)+1+x(4)-a 2*a 2*a], 'FaceColor', 'k');
rectangle('Position', [x(1)-1-x(5)-a x(2)-a 2*a 2*a], 'FaceColor', 'k');
rectangle('Position', [x(1)-a x(2)-1-x(6)-a 2*a 2*a], 'FaceColor', 'k');

% prismatic link: q2, q3, q4, q5
line([x(1)+1 x(1)+1+x(3)], [x(2) x(2)], 'LineWidth', 1);
line([x(1) x(1)], [x(2)+1 x(2)+1+x(4)], 'LineWidth', 1);
line([x(1)-1 x(1)-1-x(5)], [x(2) x(2)], 'LineWidth', 1);
line([x(1) x(1)], [x(2)-1 x(2)-1-x(6)], 'LineWidth', 1);

% display data
text(1.5, 4.5, ['q1: (', num2str(x(1), '%.2f'), ',' num2str(x(2), '%.2f'), ')']);

end