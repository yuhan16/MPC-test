function [bigMat, solverConfig] = FormulateSQP(sysdynm, sysparam)
% this function formulate the problem without substitution

% variables: X = [x, u, cn, w]
% original the model
%  J = min X'*Smat*X + 2*x0'*cmat'*X + 2*Xd*dmat'*X
%          + (xd-x0)'*Q*(xd-x0) + Xd'*Q*Xd;
%   s.t.  Amat * X = bmat*x(0),   (N*dimA x dim(X))
%         w = Dmat * x + Dlim,    (N*dimC x 1)
%         Phi(w, cn) = 0;         (N*dimC x 1)
%         xmin <= x < xim, 
%         umin <= u <= umax,  
%
% new SQP model 

% get data
A = sysdynm.A;
B = sysdynm.B;
C = sysdynm.C;
D = sysdynm.D;
dlim = sysdynm.dlim;
Q = sysdynm.Q;
Qf = sysdynm.Qf;
R = sysdynm.R;

N = sysparam.N;
xmax = sysparam.xmax;
xmin = sysparam.xmin;
umax = sysparam.umax;
umin = sysparam.umin;

dimA = size(A, 1);
dimB = size(B, 2);
dimC = size(C, 2);

% formulate big matrix
% Amat * X = bmat *x0 -> x(k+1) = Ax(k)+Bu(k)+Ccn(k)
a = kron(eye(N), -eye(dimA));
a = a + kron(tril(ones(N),-1)-tril(ones(N),-2), A);
b = kron(eye(N), B);
c = kron(eye(N), C);
bigMat.Amat = sparse([a b c zeros(N*dimA, N*dimC)]);
bigMat.bmat = sparse([-A; kron(ones(N-1, 1), zeros(dimA))]);

% w = Dmat * x + Dlim
bigMat.Dmat = sparse(kron(eye(N), D));
bigMat.Dlim = kron(ones(N, 1), dlim);

% Smat, cmat and dmat
q = blkdiag(kron(eye(N-1), Q), Qf);
r = kron(eye(N), R);
S = blkdiag(q, r);
e = 0;    % perturbed quantity so that Smat > 0
%Smat = sparse([S zeros(N*(dimA+dimB), 2*N*dimC); ...
%            zeros(2*N*dimC, N*(dimA+dimB)) zeros(2*N*dimC)]);
Smat = blkdiag(S, zeros(2*N*dimC));
bigMat.Smat = sparse( Smat + e*eye(size(Smat)) );
bigMat.cmat = zeros(N*(dimA+dimB+2*dimC), dimA);
bigMat.dmat = -[q zeros(N*dimA, N*(dimB+2*dimC))]';
bigMat.q = q;

% solver configuration 
solverConfig.senselst = repmat('=', 1, N*dimA+2*N*dimC);
solverConfig.vtypelst = repmat('C', 1, N*(dimA+dimB+2*dimC));
solverConfig.lb = [kron(ones(N,1), xmin); kron(ones(N,1), umin); zeros(2*N*dimC, 1)];
solverConfig.ub = [kron(ones(N,1), xmax); kron(ones(N,1), umax); Inf*ones(2*N*dimC, 1)];

end