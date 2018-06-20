close all
clearvars;

% sys info
dt = 0.1;
dmin = -5; 
dmax = 5;
r = 0.6;

sysparam.jointnum = 2;
sysparam.contpt = 3;
sysparam.frictdir = 2;
sysparam.xmin = [dmin; dmin; -1; -1; ];
sysparam.xmax = [dmax; dmax; 1; 1];
sysparam.umin = [-1; -1];
sysparam.umax = [1; 1];
sysparam.bigM = 200;
sysparam.iter = 100;
sysparam.N = 15;
sysparam.mu = 0.8;

% get data
jointnum = sysparam.jointnum;
contpt = sysparam.contpt;
frictdir = sysparam.frictdir;
mu = sysparam.mu;

% compute dynamics
[Jac, Msym, n, span] = ModelSymbolic1DFriction;
% mass for each link
m = [1 1];
M = [m(1) 0; 0, m(2)];
% surface normal, order: rightside, lowerside, bevelside
n = double(n);
n = n / (diag(sqrt(diag(n'*n))));   % normalize gradient 
span = span / (diag(sqrt(diag(span'*span))));

% system variables
%  - x = [x1; x2; dx1; dx2], 4x1
%  - u = [u1; u2], 2x1
%  - gamma = [cn, beta, lambda], 12x1
%  - z, bin var, 12x1

% generate system matrix
%  - x(k+1) = Ax(k) + Bu(k) + Cgamma(k)
%  - 0 <= gamma(k) _|_ Dx(k+1) + Egamma(k+1) + const >=0 
mat = [eye(jointnum), -dt*eye(jointnum); zeros(jointnum), eye(jointnum)];
sysdynm.A = inv(mat);
sysdynm.B = inv(mat) * [zeros(jointnum); inv(M)*dt];
sysdynm.C = inv(mat) * [zeros(jointnum, contpt*(2+frictdir)); ...
                        inv(M)*n inv(M)*span zeros(jointnum, contpt)];
sysdynm.D = [blkdiag(n', span'); zeros(contpt, 2*jointnum)];
%sysdynm.E = [zeros(1, 4); zeros(2, 3) ones(2, 1); sysparam.mu -ones(1,2) 0];
sysdynm.E = [zeros(contpt, contpt*(2+frictdir)); ...
         zeros(contpt*frictdir, contpt*(1+frictdir)) kron(eye(contpt), ones(frictdir, 1)); ...
         mu*eye(contpt) -kron(eye(contpt), ones(1, frictdir)) zeros(contpt)];
%sysdynm.const = [-sqrt(1)*r; zeros(3, 1)];
sysdynm.const = [-sqrt(1)*r; dmax-r; -dmin-r; zeros(contpt*(frictdir+1), 1)];
sysdynm.Q = blkdiag(10,10, 1,1);
sysdynm.Qf = sysdynm.Q;
sysdynm.R = 0.01*blkdiag(1, 1);

% get system dimension  
sysparam.dimA = size(sysdynm.A, 1);
sysparam.dimB = size(sysdynm.B, 2);
sysparam.dimC = size(sysdynm.C, 2);

% initial states and desired terminal point
sysparam.x0 = [-2; -4; 0; 0];
sysparam.xd = [4; -5; 0; 0]; 

% get data
dimA = sysparam.dimA;
dimB = sysparam.dimB;
dimC = sysparam.dimC;
N = sysparam.N;

% solve no sub problem
fprintf('\nProcessing NoSubFriction...\n');
tStart = tic; 
[bigMat_nsf, solverConfig_nsf] = FormulateNoSubFriction(sysdynm, sysparam);
tFormulate = toc(tStart); 

initcond = zeros(sysparam.N*(dimA+dimB+2*dimC), 1);
[data_nosub, datalcp] = SolveNoSubFriction(bigMat_nsf, solverConfig_nsf, sysdynm, sysparam, initcond);

fprintf('NoSubFriction done.\n');
fprintf('Formulation time for NoSub: %f s\n', tFormulate);
fprintf('Average solving time for Nosub: %f s\n', mean(data_nosub.top));

%plot_switch = [0 0 1 0 3 3];
%PlotData(datalcp, plot_switch);
plot_switch = [0 0 1 0 2 2];
PlotError({data_nosub, datalcp}, {'nosub', 'lcp'}, plot_switch)