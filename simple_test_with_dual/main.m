close all
clearvars;
% define problem 
S = blkdiag(eye(4), zeros(2));
c = [ones(4, 1); zeros(2, 1)];
N = [8/3 2; 2 5/4];
M = [2 8/3; 5/4 2];
q = [-36; -25];

epsilon = 1e-8; % stop criterion
beta = 0.7;
mu = 0.5;
delta = 1; % need to update alpha 
sigma = 0.45;    % backtracking line search param
alpha = 5;     % penalty weight coeff
rho = 0.9;      % shrinking coeff of t in backtracking line search
%initcond = [8;8; 0;0; 4/3;1];
initcond = [0;0; 0;0; q];

% split variable
x = initcond(1: 2);
y = initcond(3: 4);
w = initcond(5: 6);

% configure solver
senselst = [repmat('<', 1, 6) repmat('=', 1, 4)];
vtypelst = repmat('C', 1, 6);
%lb = zeros(6, 1);
lb = [-Inf*ones(2, 1); zeros(4, 1)];
%ub = [10*ones(4, 1); Inf*ones(2, 1)];

% formulate static matrix 
Ain = [[eye(2); -eye(2)] zeros(4); zeros(2) eye(2) zeros(2)];
bin = 10*ones(6, 1);
D = [N M -eye(2)];
H = fn_hess(S) / 2;

lambda_dual = [0;0; 0;0; 0;0];   % dual var for inequality constraint  
nu_dual = 0;
iter = 0;
step = [];
Jop = [];

while 1
    iter = iter + 1;
    if iter >= 50
        fprintf('iteration exceeds expection.\n');
        display(initcond);
        break
    end

    % formulate dynamical matrix
    f = fn_grad(initcond, c, S) + Ain'*lambda_dual;
    Phi = [phi(y(1), w(1), mu);  phi(y(2), w(2), mu)];
    [g11, g12] = phi_grad(y(1), w(1), mu);
    [g21, g22] = phi_grad(y(2), w(2), mu);
    Phi_grad = [zeros(2); blkdiag(g11, g21); blkdiag(g12, g22)];
    Aeq = [D; Phi_grad'];
    beq = [-D*initcond-q; -Phi];
    if iter == 1
        nu_dual = -inv(Aeq*Aeq') * Aeq * f;
    end
    
    try
        clear model;
        model.A = sparse([Ain; Aeq]);
        model.obj = f;
        model.Q = sparse(H);
        %model.objcon = constTerm;
        model.rhs = [bin; beq];
        model.lb = lb;
%        model.ub = ub;
        model.sense = senselst;
        model.vtype = vtypelst;
        model.modelsense = 'min';
        %model.start = initcond;
        clear params;
        params.outputflag = 0;
        
        sol = gurobi(model, params);
        
    catch gurobiError
        fprintf('Error reported\n');
        gurobiError.message
        gurobiError.message.stack
    end
    
    if ~strcmp(sol.status, 'OPTIMAL')
        fprintf('Problem %s at iter %d.\n', sol.status, iter);
        display(initcond);
        break
    end
    
    % backtracking line search
    steplen = backtrack(initcond, sol.x, c, S, D, q, mu, alpha, rho, sigma);
    step = [step steplen];
    Jop = [Jop sol.objval];
    % update penalty function
    dlambda_dual = sol.pi(1:6) - lambda_dual;
    dnu_dual = sol.pi(7: end) - nu_dual;
    if alpha <= norm(nu_dual, inf) + delta
        alpha = max([norm(nu_dual, inf)+delta alpha+2*delta]);
    end

    % update initial condition and dual variables
    initcond = initcond + steplen*sol.x;
    lambda_dual = lambda_dual + steplen*dlambda_dual;
    nu_dual = nu_dual + steplen*dnu_dual;
    mu = beta * mu;
    
    % stop creterion
    % split variable
    x = initcond(1: 2);
    y = initcond(3: 4);
    w = initcond(5: 6);
    h_stop = [D*initcond+q; phi(y(1), w(1), 0);  phi(y(2), w(2), 0)];
    dyw_stop = sol.x(3: end);
    if norm(h_stop, inf) + norm(dyw_stop, inf) <= epsilon
        fprintf('Solution is less than tolerence.\n');
        display(initcond);
        break
    end
    
end

figure;
plot(step);
figure;
plot(Jop);
