function [result, tElapsed] = FindOptimal(probParam, solverConfig, initcond)
% this function calls Gurobi and is a universal interface

% split parameters and configure Gurobi
ineqLft = probParam{1}; 
ineqRht = probParam{2};
objTerm = probParam{3};

% objective parameters
Smat = objTerm{1};
firstOrderTerm = objTerm{2};
constTerm = objTerm{3};

% matrix inequality parameters
left = []; right = [];
for i = 1: length(ineqLft)
    left = [left; ineqLft{i}];
    right = [right; ineqRht{i}];
end

% Gurobi configuration
senselst = solverConfig.senselst;
vtypelst = solverConfig.vtypelst;
lb = solverConfig.lb;
ub = solverConfig.ub;

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
    model.start = initcond;
   
    clear params;
    params.outputflag = 0;
    
    result = gurobi(model, params);   
    tElapsed = result.runtime;

catch gurobiError
    fprintf('Error reported\n');
    gurobiError.message
    gurobiError.message.stack
end

end
