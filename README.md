# MPC-test
A repository for lab code test written by Matlab. 

This repository seeks to address the robot contact optimization problem with complementarity constraints based on the idea of MPC. The problem itself can be formulated as an QPCC, and currently the big M method is proposed to solve it with a simple 1D robot model. 

Two documents illustrate the problem formulation principle and code specification.

The solver used here is Gurobi 8.0. For more details please check http://www.gurobi.com/

## Modification:

This version made a big change to the structure of the code compared with previous ones, including

1. add two formulation methods, i.e., formulation with and without substitution;

2. separate the formulation process with main test script, using cell array to convey parameters;

3. add a switch to select and test different formulation methods.

## Result:
Although the formulation method with substitution greatly reduced the number of decision variables, the solving time is nearly twice as the formulation without substitution because the matrices in the former method have higher density.  
