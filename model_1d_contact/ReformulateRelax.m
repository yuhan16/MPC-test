function [bigMatNew, solverConfigNew] = ReformulateRelax(bigMat, solverConfig, indexlst, sysparam, Dmat, Dlim)
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

% get data
dimA = sysparam.dimA;
dimB = sysparam.dimB;
dimC = sysparam.dimC;
bigM = sysparam.bigM;
N = sysparam.N;

cont = indexlst{1};
ncont = indexlst{2};
undec = indexlst{3};
contNum = length(cont);
ncontNum = length(ncont);
undecNum = length(undec);

% formulate Emat * X <= Elim & Fmat * X <= Flim, (0)
if isempty(undec)
    Emat = []; Elim = [];
    Fmat = []; Flim = [];
else
    De = zeros(undecNum, N*dimA);
    Delim = zeros(undecNum, 1);
    If = zeros(undecNum, N*dimC);
    for i = 1: undecNum
        De(i, :) = Dmat(undec(i), :);
        Delim(i) = Dlim(undec(i));
        row = zeros(1, N*dimC);
        row(undec(i)) = 1;
        If(i, :) = row;
    end
    Emat = sparse([De zeros(undecNum, N*dimB) zeros(undecNum, N*dimC) ...
                bigM*eye(undecNum)]);
    Elim = bigM*ones(undecNum, 1) - Delim;
    Fmat = sparse([zeros(undecNum, N*dimA) zeros(undecNum, N*dimB) If ...
                -bigM*eye(undecNum)]);
    Flim = zeros(undecNum, 1);
end

% formulate Gmat * X <= Glim
if isempty(contNum)
    Gmat = [];
    Glim = [];
else
    Dg = zeros(contNum, N*dimA);
    Dglim = zeros(contNum, 1);
    Ig = zeros(contNum, N*dimC);
    for i = 1: contNum
        Dg(i, :) = Dmat(cont(i), :);
        Dglim(i) = Dlim(cont(i));
        row = zeros(1, N*dimC);
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
    Dh = zeros(ncontNum, N*dimA);
    Dhlim = zeros(ncontNum, 1);
    Ih = zeros(ncontNum, N*dimC);
    for i = 1: ncontNum
        Dh(i, :) = Dmat(ncont(i), :);
        Dhlim(i) = Dlim(ncont(i));
        row = zeros(1, N*dimC);
        row(ncont(i)) = 1;
        Ih(i, :) = row;
    end
    Hmat = sparse([Dh zeros(ncontNum, N*(dimB+dimC)+undecNum); ...
                    zeros(ncontNum, N*(dimA+dimB)) Ih zeros(ncontNum, undecNum)]);
    Hlim = [bigM*ones(ncontNum, 1)-Dhlim; zeros(ncontNum, 1)];
end

% extract each term from cell array
Smat = bigMat.Smat; cmat = bigMat.cmat; dmat = bigMat.dmat;
Amat = bigMat.Amat; bmat = bigMat.bmat;
Dmat = bigMat.Dmat; Dlim = bigMat.Dlim;
% modify elements
Smat = sparse(blkdiag(Smat, zeros(undecNum)));
cmat = sparse([cmat; zeros(undecNum, dimA)]);
dmat = sparse([dmat; zeros(undecNum, N*dimA)]);
Amat = sparse([Amat zeros(N*dimA, undecNum)]);
Dmat = sparse([Dmat zeros(N*dimC, undecNum)]);
% rewrite matrix ineq
bigMatNew.Smat = Smat; bigMatNew.cmat = cmat; bigMatNew.dmat = dmat; 
bigMatNew.Amat = Amat; bigMatNew.bmat = bmat;
bigMatNew.Dmat = Dmat; bigMatNew.Dlim = Dlim;
bigMatNew.Emat = Emat; bigMatNew.Elim = Elim;
bigMatNew.Fmat = Fmat; bigMatNew.Flim = Flim; 
bigMatNew.Gmat = Gmat; bigMatNew.Glim = Glim;
bigMatNew.Hmat = Hmat; bigMatNew.Hlim = Hlim;

% solver settings
% extract cell array
senselst = solverConfig.senselst;
vtypelst = solverConfig.vtypelst;
lb = solverConfig.lb;
ub = solverConfig.ub;
% modify elements
senselst = [senselst repmat('<', 1, N*dimC)];
vtypelst = [vtypelst repmat('B', 1, undecNum)];
lb = [lb; zeros(undecNum, 1)];
ub = [ub; ones(undecNum, 1)];
% new solver configuration
solverConfigNew.senselst = senselst;
solverConfigNew.vtypelst = vtypelst;
solverConfigNew.lb = lb;
solverConfigNew.ub = ub;

end