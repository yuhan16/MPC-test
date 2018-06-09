function [data] = SolveSub(matrix_ineq, solver_setting, sysdynm, sysparam, initcond)
% this function solves the optimization problem for 'iter' times

% get data
A = sysdynm.A;
B = sysdynm.B;
C = sysdynm.C;

x0 = sysparam.x0;
xd = sysparam.xd;
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
zop = zeros(2, iter);       % the corresponding bin vars

% extract each term from cell array
Smat = matrix_ineq{1}; cmat = matrix_ineq{2}; dmat = matrix_ineq{3};
Amat = matrix_ineq{4}; 
Dmat = matrix_ineq{5}; Dlim1 = matrix_ineq{6}; Dlim2 = matrix_ineq{7};
Emat = matrix_ineq{8}; Elim1 = matrix_ineq{9}; Elim2 = matrix_ineq{10};
Fmat = matrix_ineq{11}; Flim = matrix_ineq{12};
Gmat = matrix_ineq{13}; Glim1 = matrix_ineq{14}; Glim2 = matrix_ineq{15};
Q = matrix_ineq{16}; q = matrix_ineq{17};
N = size(q, 1) / size(Q ,1);
Xd = kron(ones(N,1), xd);

% formulate Gurobi parameters
ineq_lft = {-Dmat, Emat, Fmat, Gmat};

for i = 1: iter
%    fprintf('%d iter: \n', i);
    if i == 1
        xop(:, i) = x0;
    else
        % formulate Gurobi parameters
        ineq_rht = {Dlim1+Dlim2*x0, Elim1-Elim2*x0, Flim, Glim1-Glim2*x0};
        obj_term = {Smat, 2*cmat*x0+2*dmat*Xd, (x0-xd)'*Q*(x0-xd)+(Amat*x0-Xd)'*q*(Amat*x0-Xd)};
        [result, tElapsed] = FindOptimal({ineq_lft, ineq_rht, obj_term}, solver_setting, initcond);
        X = result.x;
        % split variables
        uvar = X(1: N*dimB);
        cvar = X(N*dimB+1: N*(dimB+dimC));
        zvar = X(N*(dimB+dimC)+1: end);
        xop(:, i) = A*x0+ B*uvar(1: dimB) + C*cvar(1: dimC);
        
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
%         clf(figure(2)); 
%         PlotRobot(x0);
%         title('robot dynamics');
%         pause(0.05);
    end
    
end
data = {xop, uop, Jop, top, fop};

end