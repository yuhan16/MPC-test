function [result] = find_optimal(Amat, bmat, Cmat,  Emat, Fmat, Gmat, Smat, Q, dlim, lb, ub,...
                 bigM, senselst, vtypelst, x0, N, dimA, dimB, dimC)

xd = kron(ones(N,1),  [-2.5; 1; 0.5; 0; 0; 0]);
xd = [xd; zeros(N*(dimB+2*dimC), 1)];
try
    clear model;
    %model.A = sparse([Amat; Hmat]);
    %model.A = sparse([Amat; Hmat; Cmat]);
    model.A = sparse([Amat; -Emat; Fmat; Gmat; Cmat]);
    %model.obj = zeros(1, N*(dimA+dimB+2*dimC));
    model.obj = -2*Smat*xd;
    model.Q = Smat;
    model.objcon = x0'*Q*x0 + xd'*Smat*xd;
    %model.rhs = [bmat*x0; bigM*ones(N*dimC, 1)];
    %model.rhs = [bmat*x0; bigM*ones(N*dimC, 1); dlim];
    model.rhs = [bmat*x0; zeros(N*dimC, 1); bigM*ones(N*dimC, 1); zeros(N*dimC, 1); dlim];
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
    
    %tStart = tic; 
    result = gurobi(model, params);
    
    %tElapsed = toc(tStart) 
    %disp(result)
    %for v=1:length(names)
    %    fprintf('%s %d\n', names{v}, result.x(v));
    %end
    %fprintf('Obj: %e\n', result.objval);

catch gurobiError
    fprintf('Error reported\n');
    gurobiError.message
    gurobiError.message.stack
end

end