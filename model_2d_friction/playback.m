function playback(x, name)
% this function uses the recorded data to play back the robot dynamics
% evolution.
% input: - x: state information
%        - name: method name

videoLength = size(x, 2);
h = figure;
for i = 1: videoLength
	clf(h); 
	PlotRobot1DFriction(x(:, i));
	title(['robot dynamics (', name, ')']);
	hold on
    text(0.5, 5.5, ['Current step: ', num2str(i-1)]);
	pause(0.05);
end

%disp('Press any key to close window.');
%k = waitforbuttonpress;
%if k == 0
%    close(h);
%end

end