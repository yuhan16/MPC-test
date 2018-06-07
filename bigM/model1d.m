close all
clearvars;
% joint number
jointnum = 3;

% compute dynamics
[Jac, Msym, n] = model1dSymbolic;

% mass for each link
m = [1 1 1];
M = [m(1)+m(2)+m(3) m(2) -m(3); m(2) m(2) 0; -m(3) 0 m(3)]; 

% sys info
mu = 0.6;
dt = 0.1;
N = 20;
Q = blkdiag(10, 1, 10, 1, 1, 1);
Qf = Q;
R = blkdiag(1, 1, 1);
dmin = -5; 
dmax = 5;
xmin = [dmin; 0; 0; -1; -1; -1];
xmax = [dmax; dmax; dmax; 1; 1; 1];
umin = [-1; -1; -1];
umax = [1; 1; 1];
bigM = 10;
iter = 100;

% generate surface normal, n(:, 1) right side, n(:, 2) left side
n = double(n);

% system variables
%  - x = [x1; x2; x3; dx1; dx2; dx3], 6x1
%  - u = [u1; u2; u3], 3x1
%  - cn, 2x1
%  - z, bin var, 2x1

% generate system matrix
%  - x(k+1) = Ax(k) + Bu(k) + Ccn(k)
%  - 0 <= cn(k) _|_ Dx(k+1) + dlim >=0 
tStart = tic; 
mat = [eye(jointnum), -dt*eye(jointnum); zeros(jointnum), eye(jointnum)];
A = inv(mat);
B = inv(mat) * [zeros(jointnum); inv(M)*dt];
C = [zeros(jointnum, 2*1); inv(M)*n];
D = [n' zeros(2, jointnum)];
dlim = [dmax-1; -dmin-1];

% get system dimension  
[dimA, ~] = size(A);
[~, dimB] = size(B);
[~, dimC] = size(C);

% initial states and desired terminal point
x0 = [0; 1; 1; 0; 0; 0];
xd = [3; 1.2; 0; 0; 0; 0]; 
%xd = kron(ones(N,1),  [3; 1.2; 0; 0; 0; 0]);

ISSUB = false;   % if use method with substitution
if ISSUB==true
    [matrix_ineq, solver_setting] = FormulateSub(A, B, C, D, dlim, Q, Qf, R, N, ...
                                              bigM, xmax, xmin, umax, umin);
    initcond = zeros(N*(dimB+2*dimC), 1);
else
    [matrix_ineq, solver_setting] = FormulateNoSub(A, B, C, D, dlim, Q, Qf, R, N, ...
                                              bigM, xmax, xmin, umax, umin);
    initcond = zeros(N*(dimA+dimB+2*dimC), 1);
end
tFormulate = toc(tStart); 

% plot robot model
figure(1)
PlotRobot(x0);
title('robot model');
figure(2)

xop = zeros(dimA, iter);    % optimal states (current state)
uop = zeros(dimB, iter);    % optimal input, the first input
Jop = zeros(1, iter);       % optimal cost
top = zeros(1, iter);       % time consumed at each step
fop = zeros(2, iter);       % contact force, 1 for right, 2 for left
zop = zeros(2, iter);       % the corresponding bin vars
for i = 1: iter
    fprintf('%d iter: \n', i);
    if i == 1
        xop(:, i) = x0;
    else
        [result, tElapsed] = FindOptimal(matrix_ineq, solver_setting, x0, xd, initcond, ISSUB);
        X = result.x;
        if ISSUB==true
            % split variables
            uvar = X(1: N*dimB);
            cvar = X(N*dimB+1: N*(dimB+dimC));
            zvar = X(N*(dimB+dimC)+1: end);
            xop(:, i) = A*x0+ B*uvar(1: dimB) + C*cvar(1: dimC);
        else
            % split variables
            xvar = X(1: N*dimA+1);
            uvar = X(N*dimA+1: N*(dimA+dimB));
            cvar = X(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));
            zvar = X(N*(dimA+dimB+dimC)+1: end);
            xop(:, i) = xvar(1: dimA);
        end  
        x0 = xop(:, i);
        uop(:, i) = uvar(1: dimB);
        Jop(i) = result.objval;
        top(i) = tElapsed;
        fop(1,i) = cvar(1);
        fop(2,i) = cvar(2);
        zop(1,i) = zvar(1);
        zop(2,i) = zvar(2);
        initcond = X;
        % display solution
%         disp('current state x =');
%         disp(x0);
%         disp('optimal input u =');
%         disp(X(N*dimA+1: N*dimA+dimB));
%         disp('optimal cost J =');
%         disp(Jop(i));
    
        % plot robot dynamics
        clf(figure(2)); 
        PlotRobot(x0);
        title('robot dynamics');
        pause(0.05);
    end
    
end

% plot switch for each component: [state, input, cost, time, force]
% value = 1: plot, value = 0: not plot
plot_switch = [0 0 0 0 1];
PlotData(xop, uop, Jop, top, fop, zop, plot_switch);

% test time report
disp('%%%%%% test time report %%%%%%');
fprintf('Solve using method ');
if ISSUB
    fprintf('with substitution.\n');
else
    fprintf('without substitution. \n');
end
fprintf('problem formulation time: %f s\n', tFormulate);
fprintf('average solving time: %f s\n', mean(top));
