# MPC-test
A repository for lab code test written by Matlab.

This repository seeks to address the robot contact optimization problem with complementarity constraints based on the idea of MPC. The problem itself can be formulated as an QPCC, and currently the big M method is proposed to reformulate the problem as MIQP with a simple 1D robot model.

Two introduction document illustrates the problem formulation principle and code specification.

The solver used here is Gurobi 8.0. For more details please check http://www.gurobi.com/

## Notice
An extended 2D robot model with four possible contact and MPC without substitution are adopted in the test. Two auxiliary functions are provided: PlotData.m and playback.m.

This is a testing brunch and the more new features will be added. e.g. ploting optimal cost with respect to different initial conditions. 
