function [Jp, Jz, Jn] = findindex(y, w, mu)

Jp = [];
Jz = [];
Jn = [];

for i = 1: length(y)
    if phi(y(i), w(i), mu) > 0
        Jp = [Jp i];
    elseif phi(y(i), w(i), mu) < 0
        Jn = [Jn i];
    else
        Jz = [Jz i];
    end
end

end