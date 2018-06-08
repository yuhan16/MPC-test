function [z] = Phi(x, y, mu)

% this function formulate the big Phi function, Phi=[phi1, phi2, ...]
varlen = length(x);
z = zeros(varlen, 1);
for i = 1: varlen
    z(i) = phi(x(i), y(i), mu);
end

end