function [z] = fn(x, c, S)

z = c' * x + x' * S * x;
end