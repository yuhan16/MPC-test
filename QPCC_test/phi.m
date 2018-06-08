function [z] = phi(x, y, mu)

z = x + y - sqrt(x*x + y*y + mu);
end