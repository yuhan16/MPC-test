function [z] = theta_dot(z, varz, c, S, mu, alpha)

% split variables
x = z(1:2); 
y = z(3:4);
w = z(5:6);
dx = varz(1:2);
dy = varz(3:4);
dw = varz(5:6);
% construct index set
[Jp, Jz, Jn] = findindex(y, w, mu);
% compute each part
tmp1 = fn_grad([x; y], c, S)' * [dx; dy];
tmp2 = 0;
for i = 1: length(y)
    [phi_a, phi_b] = phi_grad(y(i), w(i), mu);
    tmp3 = alpha * [phi_a phi_b] * [dy(i); dw(i)]; 
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
        result = find(Jn == i, 1);
        if ~isempty(result)
            tmp2 = tmp2 - tmp3;
        end
    end
end
z = tmp1 + tmp2;

end