function [Jac, M, n] = model1dSymbolic
    
% joint number
jointnum = 3;

% compute Jacobian
syms x y z real
syms vx vy vz real
syms m1 m2 m3 real
m = [m1 m2 m3];
Jac = zeros(6, jointnum, jointnum);

% jacobian component: jac01, jac02, jac03
% no need for jac12 or jac13 etc.
Jac(:,:,1) = double(jacobian([0;0;vx;0;0;0], [vx;vy;vz]));
Jac(:,:,2) = double(jacobian([0;0;vx+vy;0;0;0], [vx;vy;vz]));
Jac(:,:,3) = double(jacobian([0;0;vx-vz;0;0;0], [vx;vy;vz]));

% compute Mass matrix, K = 1/2 dotq' * M * dotq
% G(q) and C(q,dotq) is 0 in this model 
M = 0;
for i = 1: jointnum
    M = M + m(i) * Jac(:,:,i)' * Jac(:,:,i);
end

Phi = [5-x-y; x-z+5];
n = [gradient(Phi(1), [x y z]) gradient(Phi(2), [x y z])];
q = [vx; vy; vz];
f = q' * M * q;
parvx = diff(f, vx);
parvy = diff(f, vy);
parvz = diff(f, vz);

end