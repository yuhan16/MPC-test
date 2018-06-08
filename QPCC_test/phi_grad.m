function [z] = phi_grad(x, y, mu)

% compute the gradient of phi(x,y)
tmp = x*x + y*x + mu;
a = 1 - x / sqrt(tmp);
b = 1 - y / sqrt(tmp);
z = [a; b];
end