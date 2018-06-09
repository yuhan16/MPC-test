function [matrix_ineq_new, solver_setting_new] = ReformulateRelax(matrix_ineq, solver_setting, indexlst, dimA, dimB, dimC, Dmat, Dlim, bigM)
% this function reformulates the relaxed problem into an MIQP so that more
% precise solution can be computed.
% new problem variables: X = [x, u, cn, zbar] (zbar undecided bin seq)
%  J = min X'*Smat*X + 2*x0'*cmat'*X + 2*Xd*dmat'*X
%          + (xd-x0)'*Q*(xd-x0) + Xd'*Q*Xd;
%   s.t.  Amat * X = bmat*x(0),   (N*dimA x dim(X))
%         Dmat * X + Dlim >= 0,   (N*dimC x dim(X))
%         Emat * X <= Elim,       (#undecided x dim(X))
%         Fmat * X <= Flim, (0)   (#undecided x dim(X))
%         Gmat * X <= Glim,       (#contact x dim(X))
%         Hmat * X <= Hlim,       (#no_contact x dim(X)) 
%         xmin <= x < xim, 
%         umin <= u <= umax, 
%         cn >=0, 

cont = indexlst{1};
ncont = indexlst{2};
undec = indexlst{3};
contNum = length(cont);
ncontNum = length(ncont);
undecNum = length(undec);
NdimA = size(Dmat, 2);
NdimC = size(Dlim, 1);
N = NdimA / dimA;

% formulate Emat * X <= Elim & Fmat * X <= Flim, (0)
if isempty(undec)
    Emat = []; Elim = [];
    Fmat = []; Flim = [];
else
    De = zeros(undecNum, NdimA);
    Delim = zeros(undecNum, 1);
    If = zeros(undecNum, NdimC);
    for i = 1: undecNum
        De(i, :) = Dmat(undec(i), :);
        Delim(i) = Dlim(undec(i));
        row = zeros(1, NdimC);
        row(undec(i)) = 1;
        If(i, :) = row;
    end
    Emat = sparse([De zeros(undecNum, N*dimB) zeros(undecNum, NdimC) ...
                bigM*eye(undecNum)]);
    Elim = bigM*ones(undecNum, 1) - Delim;
    Fmat = sparse([zeros(undecNum, NdimA) zeros(undecNum, N*dimB) If ...
                -bigM*eye(undecNum)]);
    Flim = zeros(undecNum, 1);
end

% formulate Gmat * X <= Glim
if isempty(contNum)
    Gmat = [];
    Glim = [];
else
    Dg = zeros(contNum, NdimA);
    Dglim = zeros(contNum, 1);
    Ig = zeros(contNum, NdimC);
    for i = 1: contNum
        Dg(i, :) = Dmat(cont(i), :);
        Dglim(i) = Dlim(cont(i));
        row = zeros(1, NdimC);
        row(cont(i)) = 1;
        Ig(i, :) = row;
    end
    Gmat = sparse([Dg zeros(contNum, N*(dimB+dimC)+undecNum); ...
                    zeros(contNum, N*(dimA+dimB)) Ig zeros(contNum, undecNum)]);
    Glim = [-Dglim; bigM*ones(contNum, 1)];
end

% formulate Hmat * X <= Hlim
if isempty(ncontNum)
    Hmat = [];
    Hlim = [];
else
    Dh = zeros(ncontNum, NdimA);
    Dhlim = zeros(ncontNum, 1);
    Ih = zeros(ncontNum, NdimC);
    for i = 1: ncontNum
        Dh(i, :) = Dmat(ncont(i), :);
        Dhlim(i) = Dlim(ncont(i));
        row = zeros(1, NdimC);
        row(ncont(i)) = 1;
        Ih(i, :) = row;
    end
    Hmat = sparse([Dh zeros(ncontNum, N*(dimB+dimC)+undecNum); ...
                    zeros(ncontNum, N*(dimA+dimB)) Ih zeros(ncontNum, undecNum)]);
    Hlim = [bigM*ones(ncontNum, 1)-Dhlim; zeros(ncontNum, 1)];
end

% extract each term from cell array
Smat = matrix_ineq{1}; cmat = matrix_ineq{2}; dmat = matrix_ineq{3};
Amat = matrix_ineq{4}; bmat = matrix_ineq{5};
Dmat = matrix_ineq{6}; Dlim = matrix_ineq{7};
% modify elements
Smat = sparse(blkdiag(Smat, zeros(undecNum)));
cmat = sparse([cmat; zeros(undecNum, dimA)]);
dmat = sparse([dmat; zeros(undecNum, N*dimA)]);
Amat = sparse([Amat zeros(N*dimA, undecNum)]);
Dmat = sparse([Dmat zeros(N*dimC, undecNum)]);
% rewrite matrix ineq
matrix_ineq_new = cell(1, 15);
matrix_ineq_new{1} = Smat; matrix_ineq_new{2} = cmat; matrix_ineq_new{3} = dmat; 
matrix_ineq_new{4} = Amat; matrix_ineq_new{5} = bmat;
matrix_ineq_new{6} = Dmat; matrix_ineq_new{7} = Dlim;
matrix_ineq_new{8} = Emat; matrix_ineq_new{9} = Elim;
matrix_ineq_new{10} = Fmat; matrix_ineq_new{11} = Flim; 
matrix_ineq_new{12} = Gmat; matrix_ineq_new{13} = Glim;
matrix_ineq_new{14} = Hmat; matrix_ineq_new{15} = Hlim;

% solver settings
% extract cell array
senselst = solver_setting{1};
vtypelst = solver_setting{2};
lb = solver_setting{3};
ub = solver_setting{4};
% modify elements
senselst = [senselst repmat('<', 1, N*dimC)];
vtypelst = [vtypelst repmat('B', 1, undecNum)];
lb = [lb; zeros(undecNum, 1)];
ub = [ub; ones(undecNum, 1)];
% rewrite cell array
solver_setting_new = cell(1, 4);
solver_setting_new{1} = senselst;
solver_setting_new{2} = vtypelst;
solver_setting_new{3} = lb;
solver_setting_new{4} = ub;

end