# MPC-test
A repository for lab code test written by Matlab.

This repository seeks to address the robot contact optimization problem with complementarity constraints based on the idea of MPC. The problem itself can be formulated as an QPCC, and currently the big M method is proposed to reformulate the problem as MIQP with a simple 1D robot model.

Two introduction document illustrates the problem formulation principle and code specification.

The solver used here is Gurobi 8.0. For more details please check http://www.gurobi.com/

## Notice
Three different formulation are proposed in the code, and the relaxed version is only valid for non-contact scenarios and is a good approximation for MIQP. But it will be infeasible when contact happens, which is because the relaxation method is not well considered. Anyway, it is an attempt and is kept in the code segement. For precise solution, we should run the other two formulations.

This is a developing brunch and the more new features will be added. The feature work will include adding new models, optimizing code structure and supplementing documents.
