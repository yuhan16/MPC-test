function [Jp, Jz, Jm] = FindIndex(y, w, mu)

Jp = [];
Jz = [];
Jm = [];

for i = 1: length(y)
    if phi(y(i), w(i), mu) > 0
        Jp = [Jp i];
    elseif phi(y(i), w(i), mu) < 0
        Jm = [Jm i];
    else
        Jz = [Jz i];
    end
end

end