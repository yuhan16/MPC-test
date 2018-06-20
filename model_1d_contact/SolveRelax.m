function [data] = SolveRelax(bigMat, solverConfig, sysdynm, sysparam, initcond)
% this function solves the optimization problem for 'iter' times

% get data
D = sysdynm.D;
dlim = sysdynm.dlim;

dimA = sysparam.dimA; 
dimB = sysparam.dimB;
dimC = sysparam.dimC;
x0 = sysparam.x0;
xd = sysparam.xd;
iter = sysparam.iter;
N = sysparam.N;

Smat = bigMat.Smat; cmat = bigMat.cmat; dmat = bigMat.dmat;
Amat = bigMat.Amat; bmat = bigMat.bmat;
Dmat = bigMat.Dmat; Dlim = bigMat.Dlim;
Emat = bigMat.Emat; Elim = bigMat.Elim;
Q = bigMat.Q; q = bigMat.q;
Xd = kron(ones(N,1), xd);

% allocate variables
xop = zeros(dimA, iter);    % optimal states (current state)
uop = zeros(dimB, iter);    % optimal input, the first input
Jop = zeros(1, iter);       % optimal cost
top = zeros(1, iter);       % time consumed at each step
fop = zeros(2, iter);       % contact force, 1 for right, 2 for left
%zop = zeros(2, iter);       % the corresponding bin vars

% formulate Gurobi parameters
ineqLftRlx = {Amat, -Dmat, Emat};
%figure(2)
for i = 1: iter
%    fprintf('%d iter: \n', i);
    if i == 1
        xop(:, i) = x0;
    else
        % formulate Gurobi parameters
        ineqRhtRlx = {bmat*x0, Dlim, Elim};
        objTermRlx = {Smat, 2*cmat*x0+2*dmat*Xd, (x0-xd)'*Q*(x0-xd)+Xd'*q*Xd};
        [result, tElapsed] = FindOptimal({ineqLftRlx, ineqRhtRlx, objTermRlx}, solverConfig, initcond);
        X = result.x;
        % split variables
        xvar = X(1: N*dimA);
%        uvar = X(N*dimA+1: N*(dimA+dimB));
        cvar = X(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));

        % check the solution distribution
%        check_distribution(xvar, cvar, D, dlim);
        [cont, ncont, undec] = FindContactIndex(xvar, cvar, D, dlim);
        indexlst = {cont, ncont, undec};
        [bigMatNew, solverConfigNew] = ReformulateRelax(bigMat, solverConfig, indexlst, sysparam, kron(eye(N), D), Dlim);
        % reconfigure Gurobi    
        Sm = bigMatNew.Smat; cm = bigMatNew.cmat; dm = bigMatNew.dmat;
        Am = bigMatNew.Amat; bm = bigMatNew.bmat;
        Dm = bigMatNew.Dmat; Dl = bigMatNew.Dlim;
        Em = bigMatNew.Emat; El = bigMatNew.Elim;
        Fm = bigMatNew.Fmat; Fl = bigMatNew.Flim;
        Gm = bigMatNew.Gmat; Gl = bigMatNew.Glim;
        Hm = bigMatNew.Hmat; Hl = bigMatNew.Hlim;
        
        ineqLftRef = {Am, Dm, Em, Fm, Gm, Hm};
        ineqRhtRef = {bm*x0, Dl, El, Fl, Gl, Hl};
        objTermRref = {Sm, 2*cm*x0+2*dm*Xd, (x0-xd)'*Q*(x0-xd)+Xd'*q*Xd};
        initcond = [X; 1; zeros(length(undec)-1, 1)];
        [resultRef, tElapsedRef] = FindOptimal({ineqLftRef, ineqRhtRef, objTermRref}, solverConfigNew, initcond);
        X = resultRef.x;
        % split variables
        xvar = X(1: N*dimA);
        uvar = X(N*dimA+1: N*(dimA+dimB));
        cvar = X(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));
%       zvar = X(N*(dimA+dimB+dimC)+1: end);
        
        initcond = [xvar; uvar; cvar];
        xop(:, i) = xvar(1: dimA);
        x0 = xop(:, i);
        uop(:, i) = uvar(1: dimB);
        Jop(i) = result.objval;
        top(i) = tElapsed + tElapsedRef;
        fop(1,i) = cvar(1);
        fop(2,i) = cvar(2);
%        zop(1,i) = zvar(1);
%        zop(2,i) = zvar(2);

        % plot robot dynamics
%         clf(figure(2)); 
%         PlotRobot(x0);
%         title('robot dynamics');
%         pause(0.1);
    end    
end
data = {xop, uop, Jop, top, fop};

end