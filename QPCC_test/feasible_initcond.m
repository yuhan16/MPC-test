function [initcond] = feasible_initcond(sysdynm, sysparam, bigMat)
% this function generates a strict feasible initial condition

dimA = size(sysdynm.A, 1);
dimB = size(sysdynm.B, 2);
dimC = size(sysdynm.C, 2);
A = sysdynm.A;
x0 = sysparam.x0;
N = sysparam.N;
Dmat = bigMat.Dmat;
Dlim = bigMat.Dlim;

x = zeros(N*dimA, 1);
for i = 1: N
    if i == 1
        x(1: dimA) = A*x0;
    else
        x((i-1)*dimA+1: i*dimA) = A * x((i-2)*dimA+1: (i-1)*dimA);
    end
end
w = Dmat*x + Dlim;
initcond = [x; zeros(N*(dimB+dimC), 1); w];

end