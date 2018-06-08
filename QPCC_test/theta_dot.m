function [val] = theta_dot(z, dz, bigMat, sysdynm, sysparam)
% this function computes the derivative of theta()

% get data
dimA = size(sysdynm.A, 1);
dimB = size(sysdynm.B, 2);
dimC = size(sysdynm.C, 2);
N = sysparam.N;
mu = sysparam.mu;
alpha = sysparam.alpha;
xd = sysparam.xd;
Smat = bigMat.Smat;
dmat = bigMat.dmat;

% split variables
x = z(1: N*dimA);
c = z(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));
w = z(N*(dimA+dimB+dimC)+1: end);
dx = dz(1: N*dimA);
dc = dz(N*(dimA+dimB)+1: N*(dimA+dimB+dimC));
dw = dz(N*(dimA+dimB+dimC)+1: end);

% construct index set
[Jp, Jz, Jm] = FindIndex(c, w, mu);
% compute each part
tmp1 = 2*(Smat*z+dmat*kron(ones(N,1), xd))'*dz;
tmp2 = 0;
for i = 1: length(c)
    phi_g = phi_grad(c(i), w(i), mu);
    tmp3 = alpha * [phi_g(1) phi_g(2)] * [dc(i); dw(i)]; 
    % find which set contain ith index
    if ~isempty(Jp)
        result = find(Jp == i, 1);
        if ~isempty(result)
            tmp2 = tmp2 + tmp3;
        end
    elseif ~isempty(Jz)
        result = find(Jz == i, 1);
        if ~isempty(result)
            tmp2 = tmp2 + abs(tmp3);
        end
    else
        result = find(Jm == i, 1);
        if ~isempty(result)
            tmp2 = tmp2 - tmp3;
        end
    end
end
val = tmp1 + tmp2;

end