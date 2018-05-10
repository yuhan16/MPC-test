function [result, tElapsed] = find_optimal(Amat, bmat, Dmat, Dlim, Emat, Elim, Fmat, Smat, ...
                 Q, lb, ub, senselst, vtypelst, x0, xd, N, dimC, initcond)

try
    clear model;
    model.A = sparse([Amat; -Dmat; Emat; Fmat]);
    model.obj = -2*Smat*xd;
    model.Q = Smat;
    model.objcon = x0'*Q*x0 + xd'*Smat*xd;
    model.rhs = [bmat*x0; Dlim; Elim; zeros(N*dimC, 1)];
    model.lb = lb;
    model.ub = ub;
    model.sense = senselst;
    model.vtype = vtypelst;
    model.modelsense = 'min';
    model.start = initcond;
    
    % the following properties can be omitted
    %model.start = initcond;
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
