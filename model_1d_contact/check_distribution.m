function check_distribution(xvar, cvar, D, dlim)
% this function plots the distribution of the force and the distance
% computed by relaxed model. It is used to observe the filter result with
% the inversely proportional function intuitively.
% input: - xvar: state variable
%        - cvar: contact force
%        - D: the dynamic matrix D
%        - dlim: the dynamic offset dlim

N = length(xvar) / 6;
Dmat = kron(eye(N), D);
Dlim = kron(ones(N, 1), dlim);
state = Dmat * xvar + Dlim;
h = figure;
hold on
plot(state, cvar, 'rx');
xlabel('state');
ylabel('cn');

m = 0.01: 0.01: 10;
n = 0.01 ./ m;
plot(m, n);

maxstate = max(state(:));
maxcn = max(cvar(:));
epsilon = 0.01;
threshold = 0.3;
% use inversely propotional function to select
zvar = zeros(size(cvar));
for i = 1: length(cvar)
    if state(i)*cvar(i) <= epsilon && (state(i) >= threshold*maxstate || cvar(i) >= threshold*maxcn)
        if cvar(i) >= threshold*maxcn
            zvar(i) = 1;
        else
            zvar(i) = 0;
        end
    else
        fprintf('undecided sequence: %d \n', i);
    end
end
fprintf('\n');
pause
close(h);

end