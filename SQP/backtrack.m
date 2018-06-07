function [steplen] = backtrack(z, dz, c, S, mu, alpha, rho, sigma)

dtheta = theta_dot(z, dz, c, S, mu, alpha);
ftheta = theta(z(1:2), z(3:4), z(5:6), c, S, mu, alpha);
t = 1;
while 1
    z_new = z + t*dz; 
    f1 = theta(z_new(1:2), z_new(3:4), z_new(5:6), c, S, mu, alpha);
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