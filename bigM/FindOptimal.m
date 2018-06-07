function [result, tElapsed] = FindOptimal(matrix_ineq, solver_setting, x0, xd, initcond, ISSUB)

% extract each term from cell array
Smat = matrix_ineq{1}; cmat = matrix_ineq{2}; dmat = matrix_ineq{3};
Amat = matrix_ineq{4}; bmat = matrix_ineq{5};
Dmat = matrix_ineq{6}; Dlim1 = matrix_ineq{7}; Dlim2 = matrix_ineq{8};
Emat = matrix_ineq{9}; Elim1 = matrix_ineq{10}; Elim2 = matrix_ineq{11};
Fmat = matrix_ineq{12}; Flim = matrix_ineq{13}; 
Gmat = matrix_ineq{14}; Glim1 = matrix_ineq{15}; Glim2 = matrix_ineq{16};
Q = matrix_ineq{17}; q = matrix_ineq{18};

N = size(q, 1) / size(Q ,1);

senselst = solver_setting{1};
vtypelst = solver_setting{2};
lb = solver_setting{3};
ub = solver_setting{4};

Xd = kron(ones(N,1), xd);
const_term = (x0-xd)'*Q*(x0-xd) + Xd'*q*Xd;

try
    clear model;
%    model.A = sparse([Amat; -Dmat; Emat; Fmat; Gmat]);
    model.obj = 2*cmat*x0 + 2*dmat*Xd;
    model.Q = Smat;
    if ISSUB
        model.A = sparse([-Dmat; Emat; Fmat; Gmat]);
        model.objcon = (x0-xd)'*Q*(x0-xd)+(Amat*x0-Xd)'*q*(Amat*x0-Xd);
        model.rhs = [Dlim1+Dlim2*x0; Elim1-Elim2*x0; Flim; Glim1-Glim2*x0];
    else
        model.A = sparse([Amat; -Dmat; Emat; Fmat; Gmat]);
        model.objcon = const_term;
        model.rhs = [bmat*x0; Dlim1+Dlim2*x0; Elim1-Elim2*x0; Flim];
    end 
    model.lb = lb;
    model.ub = ub;
    model.sense = senselst;
    model.vtype = vtypelst;
    model.modelsense = 'min';
    model.start = initcond;
    
    % the following properties can be omitted
    %model.varnames = [xname uname Lname zname]';
    % write model to the file
    %gurobi_write(model, 'bigM_gurobi.lp');

    clear params;
    params.outputflag = 0;
    
    tStart = tic; 
    result = gurobi(model, params);
    
    tElapsed = toc(tStart); 
    %disp(result)
    %fprintf('Obj: %e\n', result.objval);

catch gurobiError
    fprintf('Error reported\n');
    gurobiError.message
    gurobiError.message.stack
end

end
