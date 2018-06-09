function [data] = SolveRelax(matrix_ineq, solver_setting, sysdynm, sysparam, initcond)
% this function solves the optimization problem for 'iter' times

% get data
A = sysdynm.A;
B = sysdynm.B;
C = sysdynm.C;
D = sysdynm.D;
dlim = sysdynm.dlim;

x0 = sysparam.x0;
xd = sysparam.xd;
bigM = sysparam.bigM;
iter = sysparam.iter;

[dimA, ~] = size(A);
[~, dimB] = size(B);
[~, dimC] = size(C);

% allocate variables
xop = zeros(dimA, iter);    % optimal states (current state)
uop = zeros(dimB, iter);    % optimal input, the first input
Jop = zeros(1, iter);       % optimal cost
top = zeros(1, iter);       % time consumed at each step
fop = zeros(2, iter);       % contact force, 1 for right, 2 for left
%zop = zeros(2, iter);       % the corresponding bin vars

% extract each term from cell array
Smat = matrix_ineq{1}; cmat = matrix_ineq{2}; dmat = matrix_ineq{3};
Amat = matrix_ineq{4}; bmat = matrix_ineq{5};
Dmat = matrix_ineq{6}; Dlim = matrix_ineq{7};
Emat = matrix_ineq{8}; Elim = matrix_ineq{9};
Q = matrix_ineq{10}; q = matrix_ineq{11};
N = size(q, 1) / size(Q ,1);
Xd = kron(ones(N,1), xd);

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
        [result, tElapsed] = FindOptimal({ineqLftRlx, ineqRhtRlx, objTermRlx}, solver_setting, initcond);
        X = result.x;
        % split variables
        xvar = X(1: N*dimA);
%        uvar = X(N*dimA+1: N*(dimA+dimB));
        cvar = X(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));

        % check the solution distribution
%        check_distribution(xvar, cvar, D, dlim);
        [cont, ncont, undec] = FindContactIndex(xvar, cvar, D, dlim);
        indexlst = {cont, ncont, undec};
        [matrix_ineq_new, solver_setting_new] = ReformulateRelax(matrix_ineq, solver_setting, indexlst, dimA, dimB, dimC, kron(eye(N), D), Dlim, bigM);
        % reconfigure Gurobi
        Sm = matrix_ineq_new{1}; cm = matrix_ineq_new{2}; dm = matrix_ineq_new{3};
        Am = matrix_ineq_new{4}; bm = matrix_ineq_new{5};
        Dm = matrix_ineq_new{6}; Dl = matrix_ineq_new{7};
        Em = matrix_ineq_new{8}; El = matrix_ineq_new{9};
        Fm = matrix_ineq_new{10}; Fl = matrix_ineq_new{11};
        Gm = matrix_ineq_new{12}; Gl = matrix_ineq_new{13};
        Hm = matrix_ineq_new{14}; Hl = matrix_ineq_new{15};
        
        ineqLftRef = {Am, Dm, Em, Fm, Gm, Hm};
        ineqRhtRef = {bm*x0, Dl, El, Fl, Gl, Hl};
        objTermRref = {Sm, 2*cm*x0+2*dm*Xd, (x0-xd)'*Q*(x0-xd)+Xd'*q*Xd};
        initcond = [X; 1; zeros(length(undec)-1, 1)];
        [resultRef, tElapsedRef] = FindOptimal({ineqLftRef, ineqRhtRef, objTermRref}, solver_setting_new, initcond);
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