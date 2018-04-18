# MPC-test
A repository for lab code test written by Matlab.

## Modification:

Since the simple model in the simulation actually doesn't generate friction, besides there are many other methods to approximate the friction cone like Lorentz cone, we only focus on the model with contact and the optimization speed of the model via Gurobi. The friction part is deleted in this version, including

1. approximated polygonal cone matrix D1 and D2;

2. friction force components which previously are decision variables;

3. corresponding complementarity constraints about friction force;

The optimization problem becomes

![alt text](http://www.sciweavers.org/download/Tex2Img_1524085559.jpg)
