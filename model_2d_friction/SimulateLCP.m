function [sol] = SimulateLCP(u, sysdynm, sysparam, bigMat, x0)
% this function refine the solution using the princlple of LCP beacuse the
% solution of SQP method may be physically infeasible.

% get data
A = sysdynm.A;
B = sysdynm.B;
C = sysdynm.C;
D = sysdynm.D;
E = sysdynm.E;
alpha = sysdynm.const;

N = sysparam.N;
%x0 = sysparam.x0;
xd = sysparam.xd;
dimA = sysparam.dimA;
dimB = sysparam.dimB;
dimC = sysparam.dimC;

Smat = bigMat.Smat;
dmat = bigMat.dmat;
qmat = bigMat.q;
Q = bigMat.Q;

% allocate variables
xvar = zeros(N*dimA, 1);
gvar = zeros(N*dimC, 1);
wvar = zeros(N*dimC, 1);
xtmp = x0;
Xd = kron(ones(N,1), xd);

% compute a single LCP at each prediction step
for i = 1: N
    utmp = u((i-1)*dimB+1: i*dimB);
    % solve LCP
    M = D*C + E;
    q = D*A*xtmp + D*B*utmp + alpha;
    lcp = pathlcp(M, q);
    % update and record states
    xtmp = A*xtmp + B*utmp + C*lcp;
    xvar((i-1)*dimA+1: i*dimA) = xtmp;
    gvar((i-1)*dimC+1: i*dimC) = lcp;
    wvar((i-1)*dimC+1: i*dimC) = D*xtmp + E*lcp + alpha;
end
constTerm = (x0-xd)'*Q*(x0-xd)+Xd'*qmat*Xd;
cost = [xvar; u; zeros(2*N*dimC, 1)]'*Smat*[xvar; u; zeros(2*N*dimC, 1)]...
    + 2*(dmat*Xd)'*[xvar; u; zeros(2*N*dimC, 1)] ...
    + constTerm;

% assemble the solution
sol.xvar = xvar;
sol.gvar = gvar;
sol.wvar = wvar;
sol.cost = cost;

end
