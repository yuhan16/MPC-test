close all
clearvars;
% joint number
jointnum = 6;

% compute dynamics
[Jac, Msym, n] = ModelSymbolic2D;

% mass for each link
m = [2 2 1 1 1 1];
M = [m(1)+m(3)+m(4)+m(5)+m(6) 0 m(3) 0 -m(5) 0; ...
     0, m(2)+m(3)+m(4)+m(5)+m(6) 0 m(4) 0 -m(6); ...
     m(3) 0 m(3) 0 0 0; ...
     0 m(4) 0 m(4) 0 0; ...
     -m(5) 0 0 0 m(5) 0; ...
     0 -m(6) 0 0 0 m(6)];

% sys info
dt = 0.1;
dmin = -5; 
dmax = 5;

sysdynm.Q = 10*blkdiag(10,10,10,10,10,10, 10,10,10,10,10,10);
sysdynm.Qf = sysdynm.Q;
sysdynm.R = 0.01*blkdiag(1, 1, 1, 1, 1, 1);

sysparam.xmin = [dmin; dmin; 0; 0; 0; 0; -1; -1; -1; -1; -1; -1];
sysparam.xmax = [dmax; dmax; dmax; dmax; dmax; dmax; 1; 1; 1; 1; 1; 1];
sysparam.umin = [-1; -1; -2; -2; -2; -2];
sysparam.umax = [1; 1; 2; 2; 2; 2];
sysparam.bigM = 200;
sysparam.iter = 100;
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
sysdynm.C = inv(mat) * [zeros(jointnum, 4*1); inv(M)*n];
sysdynm.D = [n' zeros(4, jointnum)];
sysdynm.dlim = [dmax-1; -dmin-1; dmax-1; -dmin-1];

% get system dimension  
sysparam.dimA = size(sysdynm.A, 1);
sysparam.dimB = size(sysdynm.B, 2);
sysparam.dimC = size(sysdynm.C, 2);

% initial states and desired terminal point
sysparam.x0 = [0;0;0.5;0;1;0; 0;0;0;0.5;0.5;0];
sysparam.xd = [3;-4;0;0;0;0; zeros(jointnum, 1)]; 

% get data
dimA = sysparam.dimA;
dimB = sysparam.dimB;
dimC = sysparam.dimC;
N = sysparam.N;

% solve problem without substitution
fprintf('\nProcessing Problem...\n');
tStart = tic; 
[bigMat, solverConfig] = FormulateNoSub(sysdynm, sysparam);
tFormulate = toc(tStart); 

initcond = zeros(sysparam.N*(dimA+dimB+2*dimC), 1);
[data] = SolveNoSub(bigMat, solverConfig, sysdynm, sysparam, initcond);

fprintf('Problem solved.\n');
fprintf('Formulation time for NoSub: %f s\n', tFormulate);
fprintf('Average solving time for Nosub: %f s\n', mean(data.top));

plotSwitch = [5 3 1 1 4];
PlotData(data, plotSwitch);
