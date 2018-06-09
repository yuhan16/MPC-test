function [matrix_ineq, solver_setting] = FormulateSub(sysdynm, sysparam)   
% this function formulate the problem without substitution
% build the model
%  J = min X'*Smat*X + 2*x0'*cmat'*X + 2*Xd*dmat'*X
%          - 2Xd'*Q*Amat*x0 + (xd-x0)'*Q*(xd-x0) + Xd'*Q*Xd;
%   s.t.  Dmat * X + Dlim1 + Dlim2*x0 >= 0,   (N*dimC x dim(X))
%         Emat * X <= Elim1 - Elim2*x0,       (N*dimC x dim(X))
%         Fmat * X <= Flim (0),               (N*dimC x dim(X))
%         Gmat * X = Glim1 - Glim2*x(0),      (2N*dimA x dim(X))
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

[dimA, ~] = size(A);
[~, dimB] = size(B);
[~, dimC] = size(C);

% compute big matrix
% system dynamics: x = Amat*x0 + Bmat*u + Cmat*cn
% Notice no dynamics constraints
Amat = zeros(dimA*N, dimA);
for i = 1: N
    Amat((i-1)*dimA+1: i*dimA, :) = mpower(A, i);
end

Diag = tril(ones(N,N));
used = zeros(N,N);
for i = 1: N
    M = Diag - tril(Diag, -i) - used;
    used = used + M;
    
    if i == 1
        elmnt1 = B;
        elmnt2 = C;
        Bmat = kron(M, elmnt1);
        Cmat = kron(M, elmnt2);
    else
        elmnt1 = A * elmnt1;
        elmnt2 = A * elmnt2;
        Bmat = Bmat + kron(M, elmnt1);
        Cmat = Cmat + kron(M, elmnt2);
    end
end

% Phi(q) = [q1+q2 <= dstmax-1; q1-q3-1 >= dstmin]
% Dmat * X + Dlim >= 0 -> Phi(q(k+1)) = Dx(k+1)+dlim >= 0
d = kron(eye(N), D);
Dmat = sparse([d*Bmat d*Cmat zeros(N*dimC, N*dimC)]);
Dlim1 = kron(ones(N, 1), dlim);
Dlim2 = d*Amat;

% Emat * X <= Elim -> Phi(q(k+1)) <= M(1-z(k)) -> Dx(k+1)+dlim <= M(1-z) 
Emat = sparse([d*Bmat d*Cmat kron(eye(N), bigM*eye(dimC))]);
Elim1 = kron(ones(N, 1), [bigM; bigM]-dlim);
Elim2 = Dlim2;

% Fmat * X <= 0 -> cn(k) <= Mz(k)
Fmat = sparse([zeros(N*dimC, N*dimB) kron(eye(N), eye(dimC)) ...
               -bigM*kron(eye(N), eye(dimC))]);
Flim = zeros(N*dimC, 1);
           
% Gmat * X <= Glim1 - Glim2*x0 -> xmin <= x <= xmax
g = kron(eye(N), [eye(dimA); -eye(dimA)]);
Gmat = sparse([g*Bmat g*Cmat zeros(N*2*dimA, N*dimC)]);
Glim1 = kron(ones(N,1), [xmax; -xmin]);
Glim2 = g*Amat;

% Smat, cmat -> J = X'*Smat*X + 2*cmat'*X
q = blkdiag(kron(eye(N-1), Q), Qf);
r = kron(eye(N), R);
S = [Bmat'*q*Bmat+r Bmat'*q*Cmat; Cmat'*q*Bmat Cmat'*q*Cmat];
Smat = sparse([S zeros(N*(dimB+dimC),N*dimC); ...
               zeros(N*dimC, N*(dimB+dimC)) zeros(N*dimC)]);

c = [Amat'*q*Bmat Amat'*q*Cmat]';
cmat = [c; zeros(N*dimC, dimA)];
d = -[q*Bmat q*Cmat]';
dmat = [d; zeros(N*dimC, N*dimA)];

% write the formulation in a compact form using cell array
% notice Amat and blim are [] in this method, but we need Amat to conpute
% the objective. this Amat is not the Amat in no substitution
matrix_ineq = cell(1, 17);
matrix_ineq{1} = Smat; matrix_ineq{2} = cmat; matrix_ineq{3} = dmat;
matrix_ineq{4} = Amat;
matrix_ineq{5} = Dmat; matrix_ineq{6} = Dlim1; matrix_ineq{7} = Dlim2;
matrix_ineq{8} = Emat; matrix_ineq{9} = Elim1; matrix_ineq{10} = Elim2;
matrix_ineq{11} = Fmat; matrix_ineq{12} = Flim;
matrix_ineq{13} = Gmat; matrix_ineq{14} = Glim1; matrix_ineq{15} = Glim2;
matrix_ineq{16} = Q; matrix_ineq{17} = q; 

% solver settings
senselst = [repmat('<', 1, 3*N*dimC) repmat('<', 1, 2*N*dimA)];
vtypelst = [repmat('C', 1, N*(dimB+dimC)) repmat('B', 1, N*dimC)];
lb = [kron(ones(N,1), umin); zeros(2*N*dimC, 1)];
ub = [kron(ones(N,1), umax); Inf*ones(N*dimC, 1); ones(N*dimC, 1)];

solver_setting = cell(1, 4);
solver_setting{1} = senselst;
solver_setting{2} = vtypelst;
solver_setting{3} = lb;
solver_setting{4} = ub;

end