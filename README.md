# MPC-test
A repository for lab code test written by Matlab.

## Modification:

This sub version adds following features based on the new model.

1. add plot function for contact force and corresponding binary variables;

2. add initial condition for the optimization problem such that initial condition at step k is the optimal solution of step k-1, changing the problem into a sequential optimization problem.

The result shows the performance improves a little bit (roughly 5ms) for solving each step.
