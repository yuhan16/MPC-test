# MPC-test
A repository for lab code test written by Matlab.

## Modification:

Since the simple model in the simulation actually doesn't generate friction, besides there are many other methods to approximate the friction cone like Lorentz cone, we only focus on the model with contact and the optimization speed of the model via Gurobi. 

1. the friction part is deleted;

2. matrix notation is slightly changed, as stated in comments;

3. add a cross symbol on the base link to better identify q1. 

After modification, the optimization problem becomes

![alt text](http://www.sciweavers.org/download/Tex2Img_1524095780.jpg)
