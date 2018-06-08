function [val] = theta(z, bigMat, sysdynm, sysparam)
% this function compute the penalty function
% theta = f(x,y) + alpha(sum|phi(y,w,mu)|)

% get data
dimA = size(sysdynm.A, 1);
dimB = size(sysdynm.B, 2);
dimC = size(sysdynm.C, 2);
alpha = sysparam.alpha;
mu = sysparam.mu;
N = sysparam.N;
x0 = sysparam.x0;
xd = sysparam.xd;
Smat = bigMat.Smat;
Q = sysdynm.Q;
q = bigMat.q;
Xd = kron(ones(N,1), xd);

% split variable
x = z(1: N*dimA);
c = z(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));
w = z(N*(dimA+dimB+dimC)+1: end);

% val1 computes FB functional part, val2 computes J(x) part
val1 = alpha * sum(abs(Phi(c, w, mu)));
val2 = z'*Smat*z - 2*Xd'*q*x + Xd'*q*Xd + (x0-xd)'*Q*(x0-xd);
val = val1 + val2;

end