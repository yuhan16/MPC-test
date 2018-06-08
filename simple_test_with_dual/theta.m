function [z] = theta(z, c, S, D, q, mu, alpha)

z1 = fn(z, c, S);
z2 = [D*z+q; phi(z(3), z(5), mu); phi(z(4), z(6), mu)];

z = z1 + alpha*norm(z2, 1);
end