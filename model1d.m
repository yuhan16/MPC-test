close all
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
distmin = -5; 
distmax = 5;
xmin = [distmin; 0; 0; -1; -1; -1];
xmax = [distmax; distmax; distmax; 1; 1; 1];
umin = [-1; -1; -1];
umax = [1; 1; 1];
x0 = [0; 1; 1; 0; 0; 0];
bigM = 1000;
iter = 100;
% desired trajectory
xd = kron(ones(N,1),  [-2.5; 1; 0.5; 0; 0; 0]);
xd = [xd; zeros(N*(dimB+2*dimC), 1)];

% generate surface normal, n(:, 1) right side, n(:, 2) left side
n = double(n);

% generate D(q), isotropic friction, 
% D(q) base, D1 rightside, D2 leftside, 2fricnum # of cone edges
fricnum = 2;
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

% get system dimension  
[dimA, ~] = size(A);
[~, dimB] = size(B);
[~, dimC] = size(C);

% compute big matrix
% Amat * X = bmat *x0 -> x(k+1) = Ax(k)+Bu(k)+Clambda(k)
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
 
% cantact constraint, Phi(q) = [q1+q2<=drht; q1-q3>=dlft]
Ctmp = [1 1 0 0 0 0; -1 0 1 0 0 0];
Ctmp = kron(eye(N), Ctmp);
Cmat = [Ctmp, zeros(2*N, N*(dimB+2*dimC))];
dlim = kron(ones(N,1), [distmax; -distmin]);

% Smat -> J = X' * Smat * X
q = blkdiag(kron(eye(N-1), Q), Qf);
r = kron(eye(N), R);
S = blkdiag(q, r);
Smat = sparse([S zeros(N*(dimA+dimB), 2*N*dimC); ...
            zeros(2*N*dimC, N*(dimA+dimB)) zeros(2*N*dimC)]);

% initial term: x0'Q*x0 
%constObjTerm = x0'*Q*x0;

% build the model
%  J = min X' * Smat * X
%   s.t.  Amat * X = bmat*x(0),
%         Emat * X >= 0,
%         Fmat * X <= M,
%         Gmat * X <= 0,
%         xmin <= x < xim, 
%         umin <= u <= umax, 
%         lambda >=0, z = {0, 1} 
senselst = [repmat('=', 1, N*dimA) repmat('<', 1, 3*N*dimC+N*2)];
vtypelst = [repmat('C', 1, N*(dimA+dimB+dimC)) repmat('C', 1, N*dimC)];
lb = [kron(ones(N,1), xmin); kron(ones(N,1), umin); zeros(2*N*dimC, 1)];
ub = [kron(ones(N,1), xmax); kron(ones(N,1), umax); Inf*ones(N*dimC, 1); ones(N*dimC, 1)];
%initcond = zeros(N*(dimA+dimB+2*dimC), 1);

% plot robot model
figure(1)
plot_robot(x0);
title('robot model');
figure(2)

xop = zeros(dimA, iter);
for i = 1: iter
    if i == 1
        xop(:, i) = x0;
    else
        fprintf('\n%d iter: \n', i);
        result = find_optimal(Amat, bmat, Cmat, Emat, Fmat, Gmat, Smat, Q, lb, ub, ...
                 senselst, vtypelst, x0, xd, dlim, bigM, N, dimC);               
        X = result.x;
        xop(:, i) = X(1: dimA);
        x0 = xop(:, i);
        
        % display solution
        disp('current state x =');
        disp(x0);
        disp('optimal input u =');
        disp(X(N*dimA+1: N*dimA+dimB));
    
        % plot robot dynamics
        clf(figure(2)); 
        plot_robot(x0);
        title('robot dynamics');
        pause(0.1);
    end
    
end
        
