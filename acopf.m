% Full AC power flow equation based OPF
%   26-04-2023
clc; clear; close all;

%% load system data
define_constants;
[baseMVA, bus, gen, branch, gencost] = loadcase('case30');

%% define constants
nb = size(bus, 1);          % number of buses
nl = size(branch, 1);       % number of lines
ng = size(gen, 1);          % number of gens

Sd = bus(:, PD) + 1i * bus(:, QD);              % bus load (MVA)
Fmax = branch(:, RATE_A);                       % branch flow limit (MVA)
Vmin = bus(:, VMIN);                            % minimum voltage magnitude (p.u.)
Vmax = bus(:, VMAX);                            % maximum voltage magnitude (p.u.)
Sgmin = gen(:, PMIN) + 1i * gen(:, QMIN);       % gen max. power output (MVA)
Sgmax = gen(:, PMAX) + 1i * gen(:, QMAX);       % gen max. power output (MVA)
C2 = gencost(:, COST);                          % gen injection cost ($/MWh^2)
C1 = gencost(:, COST + 1);                      % gen injection cost ($/MWh)
C0 = gencost(:, COST + 2);                      % gen injection cost ($/h)
ref = 1;                                        % reference bus index

% -------- unit conversion --------- %
Sd = Sd / baseMVA;              % p.u.
Fmax = Fmax / baseMVA;          % p.u.
Sgmin = Sgmin / baseMVA;        % p.u.
Sgmax = Sgmax / baseMVA;        % p.u.
C2 = C2 * (baseMVA^2);          % $/h
C1 = C1 * baseMVA;              % $/h

%% connection matrix and system admittance matrix
gbus = gen(:, GEN_BUS);
Cg = sparse(gbus, (1:ng)', 1, nb, ng);

f = branch(:, F_BUS);
t = branch(:, T_BUS);
Cf = sparse((1:nl)', f, 1, nl, nb);
Ct = sparse((1:nl)', t, 1, nl, nb);

[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% ac opf
V = sdpvar(nb, 1, 'full', 'complex');   % complex bus voltage (p.u.)
Sg = sdpvar(ng, 1, 'full', 'complex');  % generator complex power injection (p.u.)

% -------- assign initial value -------- %
mpopt = mpoption( 'out.all', 0);
results = runopf('case30', mpopt);
Vm_sol = results.bus(:, VM);                    % voltage magnitude (p.u.)
Va_sol = results.bus(:, VA) * (pi/180);         % voltage angle (rad)
[Vd_sol, Vq_sol] = pol2cart(Va_sol, Vm_sol);    % complex bus voltage in cartesian coord. (p.u.)
Pg_sol = results.gen(:, PG) / baseMVA;          % generator active power injection (p.u.)
Qg_sol = results.gen(:, QG) / baseMVA;          % generator reactive power injection (p.u.)
assign(V, Vd_sol + 1i * Vq_sol);
assign(Sg, Pg_sol + 1i * Qg_sol);

Sbus = diag(V) * conj(Ybus * V);    % complex power injection at bus (p.u.)
Sf = diag(Cf * V) * conj(Yf * V);   % complex power flow of branch at "from" end (p.u.)
St = diag(Ct * V) * conj(Yt * V);   % complex power flow of branch at "to" end (p.u.)
obj = C2' * (real(Sg).^2) + C1' * real(Sg) + ones(ng, 1)' * C0;         % linear cost function ($/h)
const = [];
const = [const, Sbus + Sd - Cg * Sg == 0];                              % ac nodal power balance eqns. (p.u.)
const = [const, real(Sf).^2 + imag(Sf).^2 <= Fmax.^2];                  % branch flow limit at "from" end (p.u.)
const = [const, real(St).^2 + imag(St).^2 <= Fmax.^2];                  % branch flow limit at "to" end (p.u.)
const = [const, imag(V(ref, 1)) == 0];                                  % reference bus angle (rad)
const = [const, Vmin.^2 <= real(V).^2 + imag(V).^2 <= Vmax.^2];         % voltage magnitude (p.u.)
const = [const, Sgmin <= Sg <= Sgmax];                                  % complex power output of generator (p.u.)
opts = sdpsettings('solver', 'ipopt', 'verbose', 1, 'usex0', 1);
diag = optimize(const, obj, opts);

% -------- post-processing -------- %
V = value(V);                   % solved complex bus voltage (p.u.)
Sg = value(Sg ) * baseMVA;      % solved generator complex power injection (MVA)
cost = value(obj);              % solved generator active power injection cost ($/h)

disp('done');