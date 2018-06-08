close all
clearvars;
% joint number
jointnum = 3;

% compute dynamics
[Jac, Msym, n] = model1d_symbolic;

% mass for each link
m = [1 1 1];
M = [m(1)+m(2)+m(3) m(2) -m(3); m(2) m(2) 0; -m(3) 0 m(3)]; 

% sys info
mu = 0.6;
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
sysparam.iter = 31;
sysparam.N = 1;
sysparam.mu = 0.6;  % Fischer-Burmeiser functional param
sysparam.epsilon = 1e-9;    % stop criterion tolerance
sysparam.rho = 0.9;     % shrink param in line search
sysparam.alpha = 10;    % weight of penalty function
sysparam.sigma = 0.3;   % coeff in line search
sysparam.beta = 0.5;    % coeff for mu in FB functional
sysparam.delta = 1;

% generate surface normal, n(:, 1) right side, n(:, 2) left side
n = double(n);

% generate system matrix
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

% solve the problem
fprintf('\nProcessing SQP...\n');
tStart = tic; 
[bigMat, solverConfig] = FormulateSQP(sysdynm, sysparam);
tFormulate = toc(tStart); 
initcond = zeros(sysparam.N*(dimA+dimB+2*dimC), 1);
%initcond = feasible_initcond(sysdynm, sysparam, bigMat);
[dataSQP] = SolveSQP(bigMat, solverConfig, sysdynm, sysparam, initcond);
fprintf('SQP done.\n');
fprintf('Formulation time for SQP: %f s\n', tFormulate);
fprintf('Average solving time for SQP: %f s\n', mean(dataSQP{4}));
