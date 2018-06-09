function [result, tElapsed] = FindOptimal(matrix_ineq, solver_setting, initcond)
% this function calls Gurobi and is a universal interface

% split parameters and configure Gurobi
ineq_lft = matrix_ineq{1}; 
ineq_rht = matrix_ineq{2};
obj_term = matrix_ineq{3};

% objective parameters
Smat = obj_term{1};
first_order_term = obj_term{2};
const_term = obj_term{3};

% matrix inequality parameters
left = []; right = [];
for i = 1: length(ineq_lft)
    left = [left; ineq_lft{i}];
    right = [right; ineq_rht{i}];
end

% Gurobi configuration
senselst = solver_setting{1};
vtypelst = solver_setting{2};
lb = solver_setting{3};
ub = solver_setting{4};

try
    clear model;
    model.A = sparse(left);
    model.obj = first_order_term;
    model.Q = Smat;
    model.objcon = const_term;
    model.rhs = right;
    model.lb = lb;
    model.ub = ub;
    model.sense = senselst;
    model.vtype = vtypelst;
    model.modelsense = 'min';
    model.start = initcond;
   
    clear params;
    params.outputflag = 0;
    
%    tStart = tic; 
    result = gurobi(model, params);   
    tElapsed = result.runtime;
%    tElapsed = toc(tStart); 

catch gurobiError
    fprintf('Error reported\n');
    gurobiError.message
    gurobiError.message.stack
end

end