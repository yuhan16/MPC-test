% joint number
jointnum = 3;

% compute dynamics
[Jac, Msym, n] = model1d_symbolic;

% mass for each link
m = [1 1 1];
M = [m(1)+m(2)+m(3), m(2), -m(3); m(2), m(2), 0; ...
     -m(3), 0, m(3)]; 

% sys info
mu = 0.6;
dt = 0.1;
N = 10;
Q = eye(2*jointnum);
Qf = Q;
R = eye(jointnum);
xmin = [-5; 0; 0; -1; -1; -1];
xmax = [5; 5; 5; 1; 1; 1];
umin = [-1; -1; -1];
umax = [1; 1; 1];
x0 = [3; 1; 1; 0; 0; 0];
bigM = 1000;
iter = 100;

% generate surface normal.
% n(:, 1) represents right side, n(:, 2) left side
%n = [0, 0, -1; 0, 0, 1]';
n = double(n);

% generate D(q), isotropic friction, D(q) base, D1 rightside, D2 leftside
fricnum = 4;
dtheta = pi / fricnum;
D = zeros(jointnum, 2*fricnum);
theta = 0;
for i = 1: fricnum
    D(:, i) = [cos(theta); sin(theta); 0];
    theta = theta + dtheta;
end
D(:, fricnum+1:end) = -D(:, 1:fricnum);
D1 = roty(-90)*rotx(45)*D;
D2 = roty(135)*D;

% system variables 
%  - x = [x1; x2; x3; dx1; dx2; dx3], 6x1
%  - u = [u1; u2; u3], 3x1
%  - gamma = [cn; beta; lambda], 2x1 + 2*(2fricnum)*1 + 2x1 

% generate system matrix
%  - x(k+1) = Ax(k) + Bu(k) + Clambda(k)
%  - lambda(k) _|_ Ex(k+1) + Flmabda(k) 
mat = [eye(jointnum), -dt*eye(jointnum); zeros(jointnum), eye(jointnum)];
A = inv(mat);
B = inv(mat) * [zeros(jointnum); inv(M)*dt];
C = [zeros(jointnum, 2*(1+2*fricnum+1)); ...
      inv(M)*n, inv(M)*[D1, D2], zeros(jointnum, 2)];
E = [n'; [D1, D2]'; zeros(2, jointnum)];
E = [zeros(2+2*2*fricnum+2, 3), E];
F = [zeros(2, 2*(1+2*fricnum+1)); ...
      zeros(2*2*fricnum, 2+2*2*fricnum), kron(eye(2), ones(2*fricnum,1)); ...
      mu*eye(2), kron(eye(2), -ones(1, 2*fricnum)), zeros(2)];
 
% compute big matrix
% Amat * X = bmat *x0 -> x(k+1) = Ax(k)+Bu(k)+Clambda(k)
[dimA, ~] = size(A);
[~, dimB] = size(B);
[~, dimC] = size(C);
a = kron(eye(N), -eye(dimA));
a = a + kron(tril(ones(N),-1)-tril(ones(N),-2), A);
b = kron(eye(N), B);
c = kron(eye(N), C);
Amat = sparse([a b c zeros(N*dimA, N*dimC)]);
% bmat
bmat = sparse([-A; kron(ones(N-1, 1), zeros(dimA))]);

% Emat * X >= 0 -> Ex(k+1) + Flambda(k) >= 0
Emat = sparse([kron(eye(N), E) zeros(N*dimC, N*dimB) ...
        kron(eye(N), F) zeros(N*dimC)]);

% Fmat * X <= M -> Ex(k+1) + Flambda(k) <= M(1-z)
Fmat = sparse([kron(eye(N), E) zeros(N*dimC, N*dimB) ...
        kron(eye(N), F) kron(eye(N), bigM*eye(dimC))]);

% Gmat * X <= 0 -> lambda(k) <= Mz(k)
Gmat = sparse([zeros(N*dimC, N*dimA) zeros(N*dimC, N*dimB) ...
         kron(eye(N), eye(dimC)) -bigM*kron(eye(N), eye(dimC))]);
 
% combine Emat, Fmat and Gmat
%Hmat = -Emat+Fmat+Gmat;

% cantact constraint, Phi(q) = [q1+q2<=drht; q1-q3>=dlft]
Ctmp = [1 1 0 0 0 0; -1 0 1 0 0 0];
Ctmp = kron(eye(N), Ctmp);
Cmat = [Ctmp, zeros(2*N, N*(dimB+2*dimC))];
dlftlim = -5; drhtlim = 5;
dlim = kron(ones(N,1), [drhtlim; -dlftlim]);


% Smat -> J = X' * Smat * X
q = blkdiag(kron(eye(N-1), Q), Qf);
r = kron(eye(N), R);
S = blkdiag(q, r);
Smat = sparse([S zeros(N*(dimA+dimB), 2*N*dimC); ...
            zeros(2*N*dimC, N*(dimA+dimB)) zeros(2*N*dimC)]);

% initial term: x0'Q*x0 
constObjTerm = x0'*Q*x0;

% build the model
%  J = min X' * Smat * X
%   s.t.  Amat * X = bmat*x(0),
%         Emat * X >= 0,
%         Fmat * X <= M,
%         Gmat * X <= 0,
%         xmin <= x < xim, 
%         umin <= u <= umax, 
%         lambda >=0, z = {0, 1} 
%senselst = [repmat('=', 1, N*dimA) repmat('<', 1, N*dimC)];
%senselst = [repmat('=', 1, N*dimA) repmat('<', 1, N*dimC+N*2)];
senselst = [repmat('=', 1, N*dimA) repmat('<', 1, 3*N*dimC+N*2)];
vtypelst = [repmat('C', 1, N*(dimA+dimB+dimC)) repmat('B', 1, N*dimC)];
lb = [kron(ones(N,1), xmin); kron(ones(N,1), umin); zeros(2*N*dimC, 1)];
ub = [kron(ones(N,1), xmax); kron(ones(N,1), umax); Inf*ones(N*dimC, 1); ones(N*dimC, 1)];
initcond = zeros(N*(dimA+dimB+2*dimC), 1);

fig1 = figure(1);
axis([-5 5 -5 5]);
rectangle('Position',[x0(1)-1 -1 2 2]);
rectangle('Position',[x0(1)+1+x0(2) -0.5 1 1]);
rectangle('Position',[x0(1)-1-x0(3)-1 -0.5 1 1]);
line([x0(1)+1 x0(1)+1+x0(2)],[0 0]);
line([x0(1)-1 x0(1)-1-x0(3)],[0 0]);
figure(2)

xop = zeros(dimA, iter);
for i = 1: iter
    if i == 1
        xop(:, i) = x0;
    else
        result = find_optimal(Amat, bmat, Cmat, Emat, Fmat, Gmat, Smat, Q, dlim, lb, ub, ...
                 bigM, senselst, vtypelst, x0, N, dimA, dimB, dimC);               
        X = result.x;
        xop(:, i) = X(1: dimA);
        x0 = xop(:, i);
    fprintf('\n%d iter: ', i);
    disp(x0);
    disp(X(N*dimA+1: N*dimA+1+dimB));
    pause(0.1);
    clf(figure(2)); 
    axis([-5 5 -5 5]);
    rectangle('Position',[xop(1,i)-1 -1 2 2]);
    rectangle('Position',[xop(1,i)+1+xop(2,i) -0.5 1 1]);
    rectangle('Position',[xop(1,i)-1-xop(3,i)-1 -0.5 1 1]);
    line([xop(1,i)+1 xop(1,i)+1+xop(2,i)],[0 0]);
    line([xop(1,i)-1 xop(1,i)-1-xop(3,i)],[0 0]);
    end
    
end
        