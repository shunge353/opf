function res = dcopf(mpc, varargin)
% DC OPF on a 3-bus System
%   20-04-2023

%% load system data
if nargin ==2
    cong = 1;
else
    cong = 0;
end

define_constants;
[baseMVA, bus, gen, branch, gencost] = loadcase(mpc);

%% define constants
nb = size(bus, 1);          % number of buses
nl = size(branch, 1);       % number of lines
ng = size(gen, 1);          % number of gens

Pd = bus(:, PD);            % bus load (MW)
Fmax = branch(:, RATE_A);   % branch flow limit (MW)
Pgmin = gen(:, PMIN);       % gen max. power output (MW)
Pgmax = gen(:, PMAX);       % gen max. power output (MW)
C1 = gencost(:, COST);      % gen injection cost ($/MWh)
ref = 1;                    % reference bus index

if cong
    Fmax = [1000; 40; 1000];
end

%% build connection matrix and system susceptance matrix
gbus = gen(:, GEN_BUS);                 % connection matrix for gen
Cg = sparse(gbus, (1:ng)', 1, nb, ng);  % element i, j is 1 if gen(j) at bus i

f = branch(:, F_BUS);               % connection matrix for line and from - to buses
t = branch(:, T_BUS);               % element k, i is 1 if branch k connects "from" bus i
i = [(1:nl)'; (1:nl)'];             % element k, j is -1 if branch k connects "to" bus j
Cft = sparse(i, [f; t], [ones(nl, 1); -ones(nl, 1)], nl, nb);

b = 1 ./ branch(:, BR_X);
Bf = spdiags(b, 0, nl, nl) * Cft;   % Bf * Va is the vector of real branch power (p.u.)
Bbus = Cft' * Bf;                   % Bbus * Va is the vector of nodal real power injection (p.u.)

%% dc opf model
Pg = sdpvar(ng, 1, 'full');     % generator real power injection (MW)
Va = sdpvar(nb, 1, 'full');     % bus voltage angle (rad)

Obj = C1' * Pg;            % linear cost function ($/h)
Const = [];
Const = [Const, Bbus * Va * baseMVA + Pd - Cg * Pg == 0];   % nodal real power balance (MW)
Const = [Const, -Fmax <= Bf * Va * baseMVA <= Fmax];        % real power branch flow limit (MW)
Const = [Const, 0 <= Va(ref) <= 0];                         % reference bus (rad)
Const = [Const, Pgmin <= Pg <= Pgmax];                      % gen real power output limit (MW)
Opts = sdpsettings('solver', 'gurobi', 'verbose', 0);
diags = optimize(Const, Obj, Opts);

res.Pg = value(Pg);         % gen real power output (MW)
res.Va = value(Va);         % bus voltage angle (rad)
res.lmp = dual(Const(1));   % lmp($/MWh)
res.cost = value(Obj);      % gen real power output cost ($/h)

end