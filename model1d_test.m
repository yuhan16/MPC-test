close all
clearvars;
% joint number
jointnum = 3;

% compute dynamics
[Jac, Msym, n] = ModelSymbolic;

% mass for each link
m = [1 1 1];
M = [m(1)+m(2)+m(3) m(2) -m(3); m(2) m(2) 0; -m(3) 0 m(3)]; 

% sys info
dt = 0.1;
dmin = -5; 
dmax = 5;

sysdynm.Q = blkdiag(10, 1, 10, 1, 1, 1);
sysdynm.Qf = sysdynm.Q;
sysdynm.R = blkdiag(1, 1, 1);

sysparam.xmin = [dmin; 0; 0; -1; -1; -1];
sysparam.xmax = [dmax; dmax; dmax; 1; 1; 1];
sysparam.umin = [-1; -1; -1];
sysparam.umax = [1; 1; 1];
sysparam.bigM = 1000;
sysparam.iter = 31;
sysparam.N = 15;

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
mat = [eye(jointnum), -dt*eye(jointnum); zeros(jointnum), eye(jointnum)];
sysdynm.A = inv(mat);
sysdynm.B = inv(mat) * [zeros(jointnum); inv(M)*dt];
sysdynm.C = [zeros(jointnum, 2*1); inv(M)*n];
sysdynm.D = [n' zeros(2, jointnum)];
sysdynm.dlim = [dmax-1; -dmin-1];

% get system dimension  
dimA = size(sysdynm.A, 1);
dimB = size(sysdynm.B, 2);
dimC = size(sysdynm.C, 2);

% initial states and desired terminal point
sysparam.x0 = [2.5; 1; 1; 0; 0; 0];
sysparam.xd = [3; 1.2; 0; 0; 0; 0]; 

% solve relax problem
% fprintf('\nProcessing relax...\n');
% tStart = tic; 
% [matrix_ineq, solver_setting] = FormulateRelax(sysdynm, sysparam);
% tFormulate = toc(tStart); 
% initcond = zeros(sysparam.N*(dimA+dimB+dimC), 1);
% [data_relax] = SolveRelax(matrix_ineq, solver_setting, sysdynm, sysparam, initcond);
% fprintf('Relax done.\n');
% fprintf('Formulation time for relax: %f s\n', tFormulate);
% fprintf('Average solving time for relax: %f s\n', mean(data_relax{4}));

% solve no sub problem
fprintf('\nProcessing NoSub...\n');
tStart = tic; 
[matrix_ineq, solver_setting] = FormulateNoSub(sysdynm, sysparam);
tFormulate = toc(tStart); 
initcond = zeros(sysparam.N*(dimA+dimB+2*dimC), 1);
[data_nosub] = SolveNoSub(matrix_ineq, solver_setting, sysdynm, sysparam, initcond);
fprintf('NoSub done.\n');
fprintf('Formulation time for NoSub: %f s\n', tFormulate);
fprintf('Average solving time for Nosub: %f s\n', mean(data_nosub{4}));

% solve sub problem
fprintf('\nProcessing Sub...\n');
tStart = tic;
[matrix_ineq, solver_setting] = FormulateSub(sysdynm, sysparam);
tFormulate = toc(tStart); 
initcond = zeros(sysparam.N*(dimB+2*dimC), 1);
[data_sub] = SolveSub(matrix_ineq, solver_setting, sysdynm, sysparam, initcond);
fprintf('Sub done.\n');
fprintf('Formulation time for Sub: %f s\n', tFormulate);
fprintf('Average solving time for Sub: %f s\n', mean(data_sub{4}));

% test time report
fprintf('\n%%%%%%%%%% test report %%%%%%%%%%\n');
plot_switch = [2 0 0 1 0];
PlotError({data_sub, data_nosub}, {'Sub', 'NoSub'}, plot_switch);
