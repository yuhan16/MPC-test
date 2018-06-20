function [Jac, M, n, span] = ModelSymbolic1DFriction
% this function generates the symbolic dynamics of 1D robot model
% no rotation is assumed here 

jointnum = 2;
syms qx qy real
syms m1 m2 real
Jac = zeros(2, jointnum, jointnum);
m = [m1 m2];

% compute jacobian
Jac(:, :, 1) = double(jacobian([qx; 0], [qx; qy]));
Jac(:, :, 2) = double(jacobian([0; qy], [qx; qy]));

% compute mass matrix
M = 0;
for i = 1: jointnum
    M = M + m(i) * Jac(:, :, i)'*Jac(:, :, i);
end

% compute boundary of admissible region
r = 0.6;
Phi = [qx-qy-sqrt(2)*r; qy+5+r; -qx+5-r];
n = [gradient(Phi(1), [qx; qy]) ...
     gradient(Phi(2), [qx; qy]) ...
     gradient(Phi(3), [qx; qy])];

% compute tangent direction of each contact surface
span = [[1 -1; 1 -1] [1 -1; 0 0] [0 0; 1 -1]];

end