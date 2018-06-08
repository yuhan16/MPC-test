function [initcond, tElapsed] = FindOptimalSQP(solverParam, solverConfig, sysdynm, sysparam, bigMat, initcond)
% this function perform SQP method, including backtracking line search

%get data
ineqLft = solverParam{1};
ineqRht = solverParam{2};
objTerm = solverParam{3};
dimA = size(sysdynm.A, 1);
dimB = size(sysdynm.B, 2);
dimC = size(sysdynm.C, 2);
x0 = sysparam.x0;
N = sysparam.N;
beta = sysparam.beta;
alpha = sysparam.alpha;
mu = sysparam.mu;
delta = sysparam.delta;

Amat = bigMat.Amat;
bmat = bigMat.bmat;
Dmat = bigMat.Dmat;
Dlim = bigMat.Dlim;

% matrix inequality parameters
left = []; right = [];
for i = 1: length(ineqLft)
    left = [left; ineqLft{i}];
    right = [right; ineqRht{i}];
end

% objective parameters
Smat = objTerm{1};
firstOrderTerm = objTerm{2};
constTerm = objTerm{3};

% Gurobi configuration
senselst = solverConfig.senselst;
vtypelst = solverConfig.vtypelst;
lb = solverConfig.lb;
ub = solverConfig.ub;

iter = 0;
%tmperr = 0;
while 1
    iter = iter + 1;
    % call Gurobi
    try
        clear model;
        model.A = sparse(left);
        model.obj = firstOrderTerm;
        model.Q = Smat;
        model.objcon = constTerm;
        model.rhs = right;
        model.lb = lb;
        model.ub = ub;
        model.sense = senselst;
        model.vtype = vtypelst;
        model.modelsense = 'min';
        %model.start = initcond;
        clear params;
        params.outputflag = 0;
    
        tStart = tic; 
        sol = gurobi(model, params);
        tElapsed = toc(tStart); 
    catch gurobiError
        fprintf('Error reported\n');
        gurobiError.message
        gurobiError.message.stack
    end
    
    % split variables
%    dx = sol.x(1: N*dimA);
    du = sol.x(N*dimA+1: N*(dimA+dimB));
    dc = sol.x(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));
    dw = sol.x(N*(dimA+dimB+dimC)+1: end);
    nu_dual = sol.pi(N*(dimA+dimC)+1: end);
    
%     % check stop criterion
%     Phi_stop = Phi(dc, dw, 0);
%     if norm(Phi_stop, inf) + norm(du, inf) <= sysparam.epsilon
%         fprintf('Current solution is less than tolerance.\n');
%         display(initcond);
%         break
%     end
    
    % backtracking line search
    steplen = BacktrackLineSearch(initcond, sol.x, bigMat, sysdynm, sysparam);
   
    % update initial condition
    initcond = initcond + steplen*sol.x;
    mu = beta * mu;
    
    % check stop criterion
    Phi_stop = Phi(dc, dw, 0);
    if norm(Phi_stop, inf) + norm(du, inf) <= sysparam.epsilon
        fprintf('Current solution is less than tolerance.\n');
        display(initcond);
        break
    end

    % update penalty function
%    fprintf('error = %f\n', norm(Phi_stop, inf) + norm(dx_stop, inf));
%     if tmperr ~= 0
%         err = tmperr - norm(Phi_stop, inf) + norm(du, inf); 
%         if abs(err) <= 0.05
%             alpha = alpha + 10;
%         end
%     end
%     tmperr = norm(Phi_stop, inf) + norm(du, inf);
    if alpha <= norm(nu_dual, inf) + delta
        alpha = max([norm(nu_dual, inf)+delta alpha+2*delta]);
    end

%     % update initial condition
%     initcond = initcond + steplen*sol;
%     mu = beta * mu;
        
    % update equalities
    x = initcond(1: N*dimA);
    c = initcond(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));
    w = initcond(N*(dimA+dimB+dimC)+1: end);
    
    %ineqRht = {bmat*x0-Amat*initcond, Dmat*xtmp+Dlim-wtmp, -Phi(ctmp, wtmp, mu)};
   
    left_tmp = [zeros(N*dimC, N*(dimA+dimB)) Phi_jac(c, w, mu)];
    %right_tmp = -Phi(c, w, mu);
    left = [ineqLft{1}; ineqLft{2}; left_tmp];
    right = [bmat*x0-Amat*initcond; Dmat*x+Dlim-w; -Phi(c, w, mu)];
    
end

end