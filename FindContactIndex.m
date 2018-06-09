function [tmpc, tmpnc, tmpun] = FindContactIndex(xvar, cvar, D, dlim)
% this function returns three groups of binary var indices which represent 
% contact, not contact, undecided

% declear index vector
cont = [];
ncont = [];
undecided = [];

% get c.c. var
N = length(xvar) / 6;
Dmat = kron(eye(N), D);
Dlim = kron(ones(N, 1), dlim);
state = Dmat * xvar + Dlim;

maxstate = max(state(:));
maxcn = max(cvar(:));
epsilon = 1e-12;
threshold = 0.2;

% use inversely propotional function to determine if contact happen
for i = 1: length(cvar)
    if state(i)*cvar(i) <= epsilon && (state(i) >= threshold*maxstate || cvar(i) >= threshold*maxcn)
        if cvar(i) >= threshold*maxcn
            cont = [cont i];
        else
            ncont = [ncont i];
        end
    else
        undecided = [undecided i];
    end
end

% a new method to select contact sequence
index = [];
conttmp = [];
nconttmp = [];
undecidedtmp = [];
tmpun = [];
for i = 1: length(cvar)
    if state(i)*cvar(i) <= epsilon
        index = [index i];
    else
        undecidedtmp = [undecidedtmp i];
        tmpun = [tmpun i];
    end
end

cntmp = cvar(index);
statetmp = state(index);
maxcntmp = max(cntmp(:));
maxstatetmp = max(statetmp(:));
thcn = 0.3;
thst = 0.3;

for i = 1: length(cntmp)
    if cntmp(i) >= thcn*maxcntmp && statetmp(i) <= thst*maxstatetmp
        conttmp = [conttmp index(i)];
    elseif cntmp(i) <= thcn*maxcntmp && statetmp(i) >= thst*maxstatetmp
        nconttmp = [nconttmp index(i)];
    else
        undecidedtmp = [undecidedtmp index(i)];
    end
end


% another method to select points, using the iversely function
tmpc = [];
tmpnc = [];

th = 1e-6;
for i = 1: length(cntmp)
    if cntmp(i) >= 1e-3 && statetmp(i) <= 1e-9
        tmpc = [tmpc index(i)];
    elseif cntmp(i) <= 1e-9 && statetmp(i) >= 1e-3
        tmpnc = [tmpnc index(i)];
    else
        tmpun = [tmpun index(i)];
    end
end

end