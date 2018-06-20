function [Jac, M, n] = ModelSymbolic2D
% this function computes the sysbolic dynamics of a 2D robot
jointnum = 6;

% compute Jacobian
syms x y z real
syms q1x q1y q2 q3 q4 q5 real
syms m1 m2 m3 m4 m5 m6 real
m = [m1 m2 m3 m4 m5 m6];
Jac = zeros(2, jointnum, jointnum);

% compute jacobian
Jac(:,:,1) = double(jacobian([q1x; 0], [q1x; q1y; q2; q3; q4; q5]));
Jac(:,:,2) = double(jacobian([0; q1y], [q1x; q1y; q2; q3; q4; q5]));
Jac(:,:,3) = double(jacobian([q1x+q2; q1y], [q1x; q1y; q2; q3; q4; q5]));
Jac(:,:,4) = double(jacobian([q1x; q1y+q3], [q1x; q1y; q2; q3; q4; q5]));
Jac(:,:,5) = double(jacobian([q1x-q4; q1y], [q1x; q1y; q2; q3; q4; q5]));
Jac(:,:,6) = double(jacobian([q1x; q1y-q5], [q1x; q1y; q2; q3; q4; q5]));

% compute Mass matrix, K = 1/2 dotq' * M * dotq
% G(q) and C(q,dotq) is 0 in this model 
M = 0;
for i = 1: jointnum
    M = M + m(i) * Jac(:,:,i)' * Jac(:,:,i);
end

Phi = [-q1x-q2+5; q1x-q4+5; -q1y-q3+5; q1y-q5+5];
n = [gradient(Phi(1), [q1x; q1y; q2; q3; q4; q5]) ...
     gradient(Phi(2), [q1x; q1y; q2; q3; q4; q5]) ...
     gradient(Phi(3), [q1x; q1y; q2; q3; q4; q5]) ...
     gradient(Phi(4), [q1x; q1y; q2; q3; q4; q5])];

% forget what this code means ...
%q = [vx; vy; vz];
%f = q' * M * q;
%parvx = diff(f, vx);
%parvy = diff(f, vy);
%parvz = diff(f, vz);

end