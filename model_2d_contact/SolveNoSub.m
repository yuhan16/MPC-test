function [data] = SolveNoSub(bigMat, solverConfig, sysdynm, sysparam, initcond)
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
xop = zeros(dimA, iter+1);    % optimal states (current state)
uop = zeros(dimB, iter);    % optimal input, the first input
Jop = zeros(1, iter);       % optimal cost
top = zeros(1, iter);       % time consumed at each step
fop = zeros(dimC, iter);       % contact force
zop = zeros(dimC, iter);       % the corresponding bin vars

% formulate Gurobi parameters
ineqLft = {Amat, -Dmat, Emat, Fmat};

for i = 1: iter+1
%    fprintf('%d iter: \n', i);
    if i == 1
        xop(:, i) = x0;
    else
        % call Gurobi to solve
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

        % store and update variables 
        x0 = xop(:, i);
        uop(:, i-1) = uvar(1: dimB);
        Jop(i-1) = result.objval;
        top(i-1) = tElapsed;
        fop(:, i-1) = cvar(1: dimC);
        zop(:, i-1) = zvar(1: dimC);
        initcond = X;
        
        % plot robot dynamics
%        clf(figure(2)); 
%        PlotRobot(x0);
%        title('robot dynamics');
%        pause(0.05);
    end   
end
data.xop = xop; data.uop = uop;
data.Jop = Jop; data.top = top;
data.fop = fop; data.zop = zop;

end