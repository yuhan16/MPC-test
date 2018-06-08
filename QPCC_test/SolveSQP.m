function [data] = SolveSQP(bigMat, solverConfig, sysdynm, sysparam, initcond)
% this function solve the problem for iter times 

% get data
x0 = sysparam.x0;
xd = sysparam.xd;
iter = sysparam.iter;
N = sysparam.N;
mu = sysparam.mu;
Q = sysdynm.Q;

dimA = size(sysdynm.A, 1);
dimB = size(sysdynm.B, 2);
dimC = size(sysdynm.C, 2);

Amat = bigMat.Amat;
bmat = bigMat.bmat;
Dmat = bigMat.Dmat;
Dlim = bigMat.Dlim;
Smat = bigMat.Smat;
dmat = bigMat.dmat;
q = bigMat.q;

% allocate variables
xop = zeros(dimA, iter);    % optimal states (current state)
uop = zeros(dimB, iter);    % optimal input, the first input
Jop = zeros(1, iter);       % optimal cost
top = zeros(1, iter);       % time consumed at each step
fop = zeros(2, iter);       % contact force, 1 for right, 2 for left
wop = zeros(2, iter);       % the corresponding bin vars

Xd = kron(ones(N,1), xd);

for i = 1: iter
    %    fprintf('%d iter: \n', i);
    if i == 1
        xop(:, i) = x0;
    else
    % formulate inequalities for Gurobi
    xtmp = initcond(1: N*dimA);
    ctmp = initcond(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));
    wtmp = initcond(N*(dimA+dimB+dimC)+1: end);
    ineqLft = {Amat, [-Dmat zeros(N*dimC, N*(dimB+dimC)) eye(N*dimC)], ...
        [zeros(N*dimC, N*(dimA+dimB)) Phi_jac(ctmp, wtmp, mu)]};
    ineqRht = {bmat*x0-Amat*initcond, Dmat*xtmp+Dlim-wtmp, -Phi(ctmp, wtmp, mu)};

%    ineqLft = {Amat, [-Dmat zeros(N*dimC, N*(dimB+dimC)) eye(N*dimC)], ...
%        [zeros(N*dimC, N*(dimA+dimB)) Phi_jac(ctmp, wtmp, mu)]};
%    ineqRht = {zeros(N*dimA, 1), zeros(N*dimC, 1), -Phi(ctmp, wtmp, mu)};
    objTerm = {Smat, 2*dmat*Xd, (x0-xd)'*Q*(x0-xd)+Xd'*q*Xd};
    % call SQP method, including backtracking linear search
    [sol, tElapsed] = FindOptimalSQP({ineqLft, ineqRht, objTerm}, solverConfig, sysdynm, sysparam, bigMat, initcond);
    X = sol;
    % split variable
    xvar = X(1: N*dimA+1);
    uvar = X(N*dimA+1: N*(dimA+dimB));
    cvar = X(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));
    wvar = X(N*(dimA+dimB+dimC)+1: end);
    xop(:, i) = xvar(1: dimA); 
        
    x0 = xop(:, i);
    sysparam.x0 = x0;
    uop(:, i) = uvar(1: dimB);
%    Jop(i) = result.objval;
    top(i) = tElapsed;
    fop(1,i) = cvar(1);
    fop(2,i) = cvar(2);
    wop(1,i) = wvar(1);
    wop(2,i) = wvar(2);
    initcond = X;
    end
end
data = {xop, uop, Jop, top, fop};

end