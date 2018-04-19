function [result] = find_optimal(Amat, bmat, Cmat, Dmat, Emat, Fmat, Smat, Q, lb, ub,...
                 senselst, vtypelst, x0, xd, dlim, bigM, N, dimC)

try
    clear model;
    model.A = sparse([Amat; Cmat; -Dmat; Emat; Fmat]);
    model.obj = -2*Smat*xd;
    model.Q = Smat;
    model.objcon = x0'*Q*x0 + xd'*Smat*xd;
    model.rhs = [bmat*x0; dlim; zeros(N*dimC, 1); bigM*ones(2*N*dimC, 1)];
    model.lb = lb;
    model.ub = ub;
    model.sense = senselst;
    model.vtype = vtypelst;
    model.modelsense = 'min';
    
    % the following properties can be omitted
    %model.start = initcond;
    %model.varnames = [xname uname Lname zname]';
    % write model to the file
    %gurobi_write(model, 'bigM_gurobi.lp');

    clear params;
    params.outputflag = 0;
    
    tStart = tic; 
    result = gurobi(model, params);
    
    tElapsed = toc(tStart) 
    %disp(result)
    %fprintf('Obj: %e\n', result.objval);

catch gurobiError
    fprintf('Error reported\n');
    gurobiError.message
    gurobiError.message.stack
end

end
