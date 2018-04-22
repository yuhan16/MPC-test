close all
% joint number
jointnum = 3;

% compute dynamics
[Jac, Msym, n] = model1d_symbolic;

% mass for each link
m = [1 1 1];
M = [m(1)+m(2)+m(3) m(2) -m(3); m(2) m(2) 0; -m(3) 0 m(3)]; 

% sys info
mu = 0.6;
dt = 0.1;
N = 10;
Q = blkdiag(10, 1, 10, 1, 1, 1);
Qf = Q;
R = blkdiag(1, 1, 1);
dmin = -5; 
dmax = 5;
xmin = [dmin; 0; 0; -1; -1; -1];
xmax = [dmax; dmax; dmax; 1; 1; 1];
umin = [-1; -1; -1];
umax = [1; 1; 1];
bigM = 1000;
iter = 200;

% generate surface normal, n(:, 1) right side, n(:, 2) left side
n = double(n);

% system variables
%  - x = [x1; x2; x3; dx1; dx2; dx3], 6x1
%  - u = [u1; u2; u3], 3x1
%  - cn, 2x1
%  - z, bin var, 2x1

% generate system matrix
%  - x(k+1) = Ax(k) + Bu(k) + Ccn(k)
%  - 0 <= cn(k) _|_ Dx(k+1) + dlim >=0 
mat = [eye(jointnum), -dt*eye(jointnum); zeros(jointnum), eye(jointnum)];
A = inv(mat);
B = inv(mat) * [zeros(jointnum); inv(M)*dt];
C = [zeros(jointnum, 2*1); inv(M)*n];
D = [n' zeros(2, jointnum)];
dlim = [dmax-1; -dmin-1];

% get system dimension  
[dimA, ~] = size(A);
[~, dimB] = size(B);
[~, dimC] = size(C);

% initial states and desired terminal point
x0 = [2; 0.5; 0; 0; 0; 0];
xd = kron(ones(N,1),  [3; 1.2; 0; 0; 0; 0]);
xd = [xd; zeros(N*(dimB+2*dimC), 1)];

% compute big matrix
% Amat * X = bmat *x0 -> x(k+1) = Ax(k)+Bu(k)+Ccn(k)
a = kron(eye(N), -eye(dimA));
a = a + kron(tril(ones(N),-1)-tril(ones(N),-2), A);
b = kron(eye(N), B);
c = kron(eye(N), C);
Amat = sparse([a b c zeros(N*dimA, N*dimC)]);

% bmat
bmat = sparse([-A; kron(ones(N-1, 1), zeros(dimA))]);

% Phi(q) = [q1+q2 <= dstmax-1; q1-q3-1 >= dstmin]
% Dmat * X + Dlim >= 0 -> Phi(q(k+1)) = Dx(k+1)+dlim >= 0
Dmat = sparse([kron(eye(N), D) zeros(N*dimC, N*(dimB+2*dimC))]);
Dlim = kron(ones(N, 1), dlim);

% Emat * X <= Elim -> Phi(q(k+1)) <= M(1-z(k)) -> Dx(k+1)+dlim <= M(1-z) 
Emat = sparse([kron(eye(N), D) zeros(N*dimC, N*(dimB+dimC)) ...
               kron(eye(N), bigM*eye(dimC))]);
Elim = kron(ones(N, 1), [bigM; bigM]-dlim);

% Fmat * X <= 0 -> cn(k) <= Mz(k)
Fmat = sparse([zeros(N*dimC, N*dimA) zeros(N*dimC, N*dimB) ...
         kron(eye(N), eye(dimC)) -bigM*kron(eye(N), eye(dimC))]);

% Smat -> J = X' * Smat * X
q = blkdiag(kron(eye(N-1), Q), Qf);
r = kron(eye(N), R);
S = blkdiag(q, r);
Smat = sparse([S zeros(N*(dimA+dimB), 2*N*dimC); ...
            zeros(2*N*dimC, N*(dimA+dimB)) zeros(2*N*dimC)]);

% initial term: x0'Q*x0 
initTerm = x0'*Q*x0;

% build the model
%  J = min X' * Smat * X + x0'*Q*x0;
%   s.t.  Amat * X = bmat*x(0),   (N*dimA x dim(X))
%         Dmat * X + Dlim >= 0,   (N*dimC x dim(X))
%         Emat * X <= Elim,       (N*dimC x dim(X))
%         Fmat * X <= 0,          (N*dimC x dim(X))
%         xmin <= x < xim, 
%         umin <= u <= umax, 
%         cn >=0, z = {0, 1} 

senselst = [repmat('=', 1, N*dimA) repmat('<', 1, 3*N*dimC)];
vtypelst = [repmat('C', 1, N*(dimA+dimB+dimC)) repmat('B', 1, N*dimC)];
lb = [kron(ones(N,1), xmin); kron(ones(N,1), umin); zeros(2*N*dimC, 1)];
ub = [kron(ones(N,1), xmax); kron(ones(N,1), umax); Inf*ones(N*dimC, 1); ones(N*dimC, 1)];
%initcond = zeros(N*(dimA+dimB+2*dimC), 1);

% plot robot model
figure(1)
plot_robot(x0);
title('robot model');
figure(2)

xop = zeros(dimA, iter);
uop = zeros(dimB, iter);
Jop = zeros(1, iter);
top = zeros(1, iter);
for i = 1: iter
    fprintf('%d iter: \n', i);
    if i == 1
        xop(:, i) = x0;
    else
        [result, tElapsed] = find_optimal(Amat, bmat, Dmat, Dlim, Emat, Elim, Fmat, Smat, ...
                  Q, lb, ub, senselst, vtypelst, x0, xd, N, dimC);

        X = result.x;
        xop(:, i) = X(1: dimA);
        x0 = xop(:, i);
        uop(:, i) = X(dimA+1: dimA+dimB);
        Jop(i) = result.objval;
        top(i) = tElapsed;
        
        % display solution
%         disp('current state x =');
%         disp(x0);
%         disp('optimal input u =');
%         disp(X(N*dimA+1: N*dimA+dimB));
%         disp('optimal cost J =');
%         disp(Jop(i));
    
        % plot robot dynamics
        clf(figure(2)); 
        plot_robot(x0);
        title('robot dynamics');
        pause(0.05);
    end
    
end

% plot switch for each component: [state, input, cost, time]
% value = 1: plot, value = 0: not plot
plot_switch = [1 0 0 1];
plot_data(xop, uop, Jop, top, plot_switch);
