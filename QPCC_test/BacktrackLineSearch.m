function steplen = BacktrackLineSearch(z, dz, bigMat, sysdynm, sysparam)
% this function implement the backtracking line search

t = 1;
rho = sysparam.rho;
sigma = sysparam.sigma;

% pre compute 2 terms
ftheta = theta(z, bigMat, sysdynm, sysparam);
dtheta = theta_dot(z, dz, bigMat, sysdynm, sysparam);
while 1
    z_new = z + t*dz;
    f1 = theta(z_new, bigMat, sysdynm, sysparam);
    f2 = sigma * t * dtheta;
    
    if f1 <= f2 + ftheta
        steplen = t;
        break
    else
        t = rho * t;
    end
end

end