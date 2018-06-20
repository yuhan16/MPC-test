# MPC-test
A repository for lab code test written by Matlab.

This repository seeks to address the robot contact optimization problem with complementarity constraints based on the idea of MPC. The problem itself can be formulated as an QPCC, and currently the big M method is proposed to reformulate the problem as MIQP with different simple robot models.

The introduction documents illustrating the problem formulation principle and code specification are included in each folder. 

The solver used here is Gurobi 8.0. For more details please check http://www.gurobi.com/

## Notice
This is a developing brunch. More new features will be added including adding new models, optimizing code structure and supplementing documents.

## Supplement
### 1D Model with frictionless contact
Three different formulation are proposed in this part. The relaxed version is only valid for non-contact scenarios and is a good approximation for MIQP, but it becomes infeasible when contact happens, which is because the relaxation method is not well considered. Anyway, it is an attempt and is kept in the code segement. For precise solution, you should run the other two formulations.

### 2D Model with frictionless contact
An extended 2D robot model with four possible contact and MPC without substitution are adopted in this part. Two auxiliary functions are provided: PlotData.m and playback.m.
