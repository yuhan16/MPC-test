function PlotRobot( x0 )

% define axis
axis([-6 6 -6 6]);
hold on
% barrier
line([-5,-5], [-6,6], 'Color','black', 'LineWidth', 1);
line([5,5], [-6,6], 'Color','black', 'LineWidth', 1);

% robot link
% base link
rectangle('Position',[x0(1)-1 -1 2 2]);
plot(x0(1), 0, 'kx');
% end effector
a = 0.15;
rectangle('Position',[x0(1)+1+x0(2)-a -a 2*a 2*a], 'FaceColor', 'k');
rectangle('Position',[x0(1)-1-x0(3)-a -a 2*a 2*a], 'FaceColor', 'k');
% prismatic link
line([x0(1)+1 x0(1)+1+x0(2)],[0 0], 'LineWidth', 1);
line([x0(1)-1 x0(1)-1-x0(3)],[0 0], 'LineWidth', 1);

end