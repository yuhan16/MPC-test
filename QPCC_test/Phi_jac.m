function [z] = Phi_jac(x, y, mu)

% this function computes the Jacobian of Phi(x,y,mu)
% the dimension of Jacobian is NdimC x 2NdimC, J = [D1, D2], where D1, D2
% are diagonal matrices

varlen = length(x);
d1 = zeros(varlen);
d2 = zeros(varlen);

for i = 1: varlen
    tmp = phi_grad(x(i), y(i), mu);
    d1(i, i) = tmp(1);
    d2(i ,i) = tmp(2);
end
z = [d1 d2];

end