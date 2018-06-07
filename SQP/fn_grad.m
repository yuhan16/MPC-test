function [y] = fn_grad(x, c, S)

y = c + 2*S*x;
%y = [y; zeros(2, 1)];
end