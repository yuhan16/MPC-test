function [phi_a, phi_b] = phi_grad(x, y, mu)

tmp = x^2 + y^2 + mu;
phi_a = 1 - x / sqrt(tmp);
phi_b = 1 - y / sqrt(tmp);
end