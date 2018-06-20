function playback(x, name)
% this function uses the recorded data to play back the robot dynamics
% evolution.
% input: - x: state information
%        - name: method name

videoLength = size(x, 2);
h = figure;
for i = 1: videoLength
	clf(h); 
	PlotRobot2D(x(:, i));
	title(['robot dynamics (', name, ')']);
	hold on
    text(1, 5.5, ['Current step: ', num2str(i)]);
	pause(0.1);
end

disp('Click window to close.');
k = waitforbuttonpress;
if k == 0
    close(h);
end

end