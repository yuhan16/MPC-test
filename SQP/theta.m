function [z] = theta(x, y, w, c, S, mu, alpha)

z1 = fn([x; y], c, S);
z2 = alpha * (abs(phi(y(1), w(1), mu)) + abs(phi(y(2), w(2), mu)));
z = z1 + z2;
end