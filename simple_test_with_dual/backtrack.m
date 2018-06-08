function [steplen] = backtrack(z, dz, c, S, D, q, mu, alpha, rho, sigma)

%dtheta = theta_dot(z, dz, c, S, mu, alpha);
dtheta = fn_grad(z, c, S)'*dz ...
         - alpha*norm([D*z+q; phi(z(3), z(5), mu); phi(z(4), z(6), mu)], 1);
ftheta = theta(z, c, S, D, q, mu, alpha);
t = 1;
while 1
    z_new = z + t*dz; 
    f1 = theta(z_new, c, S, D, q, mu, alpha);
%    f2 = theta(z(1:2), z(3:4), z(5:6), c, S, mu, alpha);
    f3 = sigma * t * dtheta;
    if f1 <= ftheta + f3
        steplen = t;
        break
    else
        t = rho * t;
    end
end

end