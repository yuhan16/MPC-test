function [bigMat, solverConfig] = FormulateNoSub(sysdynm, sysparam)
% this function formulate the problem without substitution
% build the model
%  J = min X'*Smat*X + 2*x0'*cmat'*X + 2*Xd*dmat'*X
%          + (xd-x0)'*Q*(xd-x0) + Xd'*Q*Xd;
%   s.t.  Amat * X = bmat*x(0),   (N*dimA x dim(X))
%         Dmat * X + Dlim >= 0,   (N*dimC x dim(X))
%         Emat * X <= Elim,       (N*dimC x dim(X))
%         Fmat * X <= Flim (0),   (N*dimC x dim(X))
%         xmin <= x < xim, 
%         umin <= u <= umax, 
%         cn >=0, 
%         z = {0, 1} 

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
bigM = sysparam.bigM;
xmax = sysparam.xmax;
xmin = sysparam.xmin;
umax = sysparam.umax;
umin = sysparam.umin;

dimA = sysparam.dimA;
dimB = sysparam.dimB;
dimC = sysparam.dimC;

% compute big matrix
% Amat * X = bmat *x0 -> x(k+1) = Ax(k)+Bu(k)+Ccn(k)
a = kron(eye(N), -eye(dimA));
a = a + kron(tril(ones(N),-1)-tril(ones(N),-2), A);
b = kron(eye(N), B);
c = kron(eye(N), C);
bigMat.Amat = sparse([a b c zeros(N*dimA, N*dimC)]);

% bmat
bigMat.bmat = sparse([-A; kron(ones(N-1, 1), zeros(dimA))]);

% Phi(q) = [q1+q2 <= dstmax-1; q1-q3-1 >= dstmin]
% Dmat * X + Dlim >= 0 -> Phi(q(k+1)) = Dx(k+1)+dlim >= 0
bigMat.Dmat = sparse([kron(eye(N), D) zeros(N*dimC, N*(dimB+2*dimC))]);
bigMat.Dlim = kron(ones(N, 1), dlim);

% Emat * X <= Elim -> Phi(q(k+1)) <= M(1-z(k)) -> Dx(k+1)+dlim <= M(1-z) 
bigMat.Emat = sparse([kron(eye(N), D) zeros(N*dimC, N*(dimB+dimC)) ...
               kron(eye(N), bigM*eye(dimC))]);
bigMat.Elim = kron(ones(N, 1), [bigM; bigM]-dlim);

% Fmat * X <= 0 -> cn(k) <= Mz(k)
bigMat.Fmat = sparse([zeros(N*dimC, N*dimA) zeros(N*dimC, N*dimB) ...
         kron(eye(N), eye(dimC)) -bigM*kron(eye(N), eye(dimC))]);
bigMat.Flim = zeros(N*dimC, 1);

% Smat -> J = X' * Smat * X
q = blkdiag(kron(eye(N-1), Q), Qf);
r = kron(eye(N), R);
S = blkdiag(q, r);
bigMat.Smat = sparse([S zeros(N*(dimA+dimB), 2*N*dimC); ...
            zeros(2*N*dimC, N*(dimA+dimB)) zeros(2*N*dimC)]);
bigMat.cmat = zeros(N*(dimA+dimB+2*dimC), dimA);
bigMat.dmat = -[q zeros(N*dimA, N*(dimB+2*dimC))]';
bigMat.Q = Q;
bigMat.q = q;

% solver settings
solverConfig.senselst = [repmat('=', 1, N*dimA) repmat('<', 1, 3*N*dimC)];
solverConfig.vtypelst = [repmat('C', 1, N*(dimA+dimB+dimC)) repmat('B', 1, N*dimC)];
solverConfig.lb = [kron(ones(N,1), xmin); kron(ones(N,1), umin); zeros(2*N*dimC, 1)];
solverConfig.ub = [kron(ones(N,1), xmax); kron(ones(N,1), umax); Inf*ones(N*dimC, 1); ones(N*dimC, 1)];

end