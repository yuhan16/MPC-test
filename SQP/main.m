clearvars;
% define problem 
S = eye(4);
c = ones(4, 1);
N = [8/3 2; 2 5/4];
M = [2 8/3; 5/4 2];
q = [-36; -25];

% configure solver
senselst = repmat('=', 1, 4);
vtypelst = repmat('C', 1, 6);
lb = zeros(6, 1);
ub = [10*ones(4, 1); Inf*ones(2, 1)];

epsilon = 1e-9; % stop criterion
beta = 0.7; 
mu = 0.5;
delta = 1; % need to update alpha 
sigma = 0.3;    % backtracking line search param
alpha = 5;     % penalty weight coeff
rho = 0.9;      % shrinking coeff of t in backtracking line search
%initcond = [8;8; 0;0; 4/3;1];
initcond = [0;0; 0;0; q];
iter = 0;
step = [];
Jop = [];

% split variable
x = initcond(1: 2);
y = initcond(3: 4);
w = initcond(5: 6);

while 1
    iter = iter + 1;
    if iter >= 50
        fprintf('iteration exceeds expection.\n');
        display(initcond);
        break
    end
%     % split variable
%     x = initcond(1: 2);
%     y = initcond(3: 4);
%     w = initcond(5: 6);
    % formulate the problem using quadprog()
    f = [fn_grad([x; y], c, S); zeros(2, 1)];
    H = fn_hess(S);
    Phi = [phi(y(1), w(1), mu);  phi(y(2), w(2), mu)];
    [g11, g12] = phi_grad(y(1), w(1), mu);
    [g21, g22] = phi_grad(y(2), w(2), mu);
    Phi_grad = [blkdiag(g11, g21); blkdiag(g12, g22)];
    Aeq = [N M -eye(2); zeros(2) Phi_grad'];
    %beq = [zeros(2,1); -Phi];
    beq = [w-N*x-M*y-q; -Phi];
    
    try
        clear model;
        model.A = sparse(Aeq);
        model.obj = f;
        model.Q = sparse(H);
        %model.objcon = constTerm;
        model.rhs = beq;
        model.lb = lb;
        model.ub = ub;
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
    steplen = backtrack(initcond, sol.x, c, S, mu, alpha, rho, sigma);
    step = [step steplen];
    Jop = [Jop sol.objval];
    % update penalty function
    nu_dual = sol.pi(3: 4);
    if alpha <= norm(nu_dual, inf) + delta
        alpha = max([norm(nu_dual, inf)+delta alpha+2*delta]);
    end

    % update initial condition
    initcond = initcond + steplen*sol.x;
    mu = beta * mu;
    
    % stop creterion
    % split variable
    x = initcond(1: 2);
    y = initcond(3: 4);
    w = initcond(5: 6);
    Phi_stop = [phi(y(1), w(1), 0);  phi(y(2), w(2), 0)];
    dx_stop = sol.x(1:2);
    if norm(Phi_stop, inf) + norm(dx_stop, inf) <= epsilon
        fprintf('Solution is less than tolerence.\n');
        display(initcond);
        break
    end
    
end

figure;
plot(step);
figure;
plot(Jop);
