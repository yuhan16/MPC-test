function [matrix_ineq, solver_setting] = FormulateRelax(sysdynm, sysparam)
% this function relaxs the c.c. using linear inequalities
% build the model, variable X = [x, u, cn]
%  J = min X'*Smat*X + 2*x0'*cmat'*X + 2*Xd*dmat'*X
%          + (xd-x0)'*Q*(xd-x0) + Xd'*Q*Xd;
%   s.t.  Amat * X = bmat*x(0),   (N*dimA x dim(X))
%         Dmat * X + Dlim >= 0,   (N*dimC x dim(X))
%         Emat * X <= Elim,       (N*dimC x dim(X))
%         xmin <= x < xim, 
%         umin <= u <= umax, 
%         cn >=0, 

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

[dimA, ~] = size(A);
[~, dimB] = size(B);
[~, dimC] = size(C);

% compute big matrix
% Amat * X = bmat *x0 -> x(k+1) = Ax(k)+Bu(k)+Ccn(k)
a = kron(eye(N), -eye(dimA));
a = a + kron(tril(ones(N),-1)-tril(ones(N),-2), A);
b = kron(eye(N), B);
c = kron(eye(N), C);
Amat = sparse([a b c]);

% bmat
bmat = sparse([-A; kron(ones(N-1, 1), zeros(dimA))]);

% Phi(q) = [q1+q2 <= dstmax-1; q1-q3-1 >= dstmin]
% Dmat * X + Dlim >= 0 -> Phi(q(k+1)) = Dx(k+1)+dlim >= 0
Dmat = sparse([kron(eye(N), D) zeros(N*dimC, N*(dimB+dimC))]);
Dlim = kron(ones(N, 1), dlim);

% Emat * X <= Elim -> Phi(q(k+1)) <= M(1-z(k)) -> Dx(k+1)+dlim <= M(1-z) 
Emat = sparse([kron(eye(N), D) zeros(N*dimC, N*dimB) ...
               kron(eye(N), bigM*eye(dimC))]);
Elim = kron(ones(N, 1), [bigM; bigM]) - Dlim;

% Smat -> J = X' * Smat * X
q = blkdiag(kron(eye(N-1), Q), Qf);
r = kron(eye(N), R);
S = blkdiag(q, r);
Smat = sparse([S zeros(N*(dimA+dimB), N*dimC); ...
            zeros(N*dimC, N*(dimA+dimB)) zeros(N*dimC)]);
cmat = sparse(zeros(N*(dimA+dimB+dimC), dimA));
dmat = sparse(-[q zeros(N*dimA, N*(dimB+dimC))]');

% write the formulation in a compact form using cell array
matrix_ineq = cell(1, 11);
matrix_ineq{1} = Smat; matrix_ineq{2} = cmat; %zeros(N*(dimA+dimB+2*dimC), dimA);
matrix_ineq{3} = dmat; %-[q zeros(N*dimA, N*(dimB+2*dimC))]';
matrix_ineq{4} = Amat; matrix_ineq{5} = bmat;
matrix_ineq{6} = Dmat; matrix_ineq{7} = Dlim;
matrix_ineq{8} = Emat; matrix_ineq{9} = Elim;
matrix_ineq{10} = Q; matrix_ineq{11} = q;
        
% solver settings
senselst = [repmat('=', 1, N*dimA) repmat('<', 1, 2*N*dimC)];
vtypelst = repmat('C', 1, N*(dimA+dimB+dimC));
lb = [kron(ones(N,1), xmin); kron(ones(N,1), umin); zeros(N*dimC, 1)];
ub = [kron(ones(N,1), xmax); kron(ones(N,1), umax); Inf*ones(N*dimC, 1)];

solver_setting = cell(1, 4);
solver_setting{1} = senselst;
solver_setting{2} = vtypelst;
solver_setting{3} = lb;
solver_setting{4} = ub;

end