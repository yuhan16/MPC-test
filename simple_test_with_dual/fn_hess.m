function [H] = fn_hess(S)

e = 0;
H = 2 * S;
H = H + e*eye(size(H));
end