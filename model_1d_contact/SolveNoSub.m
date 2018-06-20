function [data] = SolveNoSub(bigMat, solverConfig, sysparam, initcond)
% this function solves the optimization problem for 'iter' times

% get data
dimA = sysparam.dimA;
dimB = sysparam.dimB;
dimC = sysparam.dimC;

x0 = sysparam.x0;
xd = sysparam.xd;
iter = sysparam.iter;
N = sysparam.N;

Smat = bigMat.Smat; dmat = bigMat.dmat;
Amat = bigMat.Amat; bmat = bigMat.bmat;
Dmat = bigMat.Dmat; Dlim = bigMat.Dlim;
Emat = bigMat.Emat; Elim = bigMat.Elim;
Fmat = bigMat.Fmat; Flim = bigMat.Flim;
Q = bigMat.Q; q = bigMat.q;
Xd = kron(ones(N,1), xd);

% allocate variables
xop = zeros(dimA, iter);    % optimal states (current state)
uop = zeros(dimB, iter);    % optimal input, the first input
Jop = zeros(1, iter);       % optimal cost
top = zeros(1, iter);       % time consumed at each step
fop = zeros(2, iter);       % contact force, 1 for right, 2 for left
zop = zeros(2, iter);       % the corresponding bin vars

% formulate Gurobi parameters
ineqLft = {Amat, -Dmat, Emat, Fmat};

for i = 1: iter
%    fprintf('%d iter: \n', i);
    if i == 1
        xop(:, i) = x0;
    else
        % formulate Gurobi parameters
        ineqRht = {bmat*x0, Dlim, Elim, Flim};
        objTerm = {Smat, 2*dmat*Xd, (x0-xd)'*Q*(x0-xd)+Xd'*q*Xd};
        [result, tElapsed] = FindOptimal({ineqLft, ineqRht, objTerm}, solverConfig, initcond);
        X = result.x;
        % split variable
        xvar = X(1: N*dimA+1);
        uvar = X(N*dimA+1: N*(dimA+dimB));
        cvar = X(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));
        zvar = X(N*(dimA+dimB+dimC)+1: end);
        xop(:, i) = xvar(1: dimA); 
        
        x0 = xop(:, i);
        uop(:, i) = uvar(1: dimB);
        Jop(i) = result.objval;
        top(i) = tElapsed;
        fop(1,i) = cvar(1);
        fop(2,i) = cvar(2);
        zop(1,i) = zvar(1);
        zop(2,i) = zvar(2);
        initcond = X;

        % plot robot dynamics
%        clf(figure(2)); 
%        PlotRobot(x0);
%        title('robot dynamics');
%        pause(0.05);
    end   
end
data = {xop, uop, Jop, top, fop};

end