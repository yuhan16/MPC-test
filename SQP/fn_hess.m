function [H] = fn_hess(S)

e = 0;
H = blkdiag(2*S, zeros(2));
H = H + e*eye(size(H));
end