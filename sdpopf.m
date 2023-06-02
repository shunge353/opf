% Examing the Limits of the Application of Semidefinite Programming
% to Power Flow Problems
%
%   DOI: 10.1109/Allerton.2011.6120344
%
clc; clear; close all;

%% load system data
define_constants;
[baseMVA, bus, gen, branch, gencost] = loadcase('pglib_opf_case3_lmbd');

branch(:, branch(:,BR_R)==0) = 1e-5;  % make the graph of Re(Y) connected

nlmp = true;    % make line flow limit stressed
% nlmp = false;   % make line flow limit stressed
if ~nlmp
    branch(2, RATE_A) = 60;
end

nb = size(bus,1);       % number of buses
nl = size(branch,1);    % number of branches
ng = size(gen,1);       % number of gens

C2 = gencost(:, COST);      % gen quadratic cost parameter ($/MWh^2)
C1 = gencost(:, COST+1);    % gen linear cost parameter ($/MWh)
C0 = gencost(:, COST+2);    % gen constant cost parameter ($/h)

gbus = gen(:,GEN_BUS);
Cg = sparse(gbus, (1:ng)', 1, nb, ng);

Pgmin = gen(:, PMIN);   % gen min active power output (MW)
Pgmax = gen(:, PMAX);   % gen max active power output (MW)
Qgmin = gen(:, QMIN);   % gen min reactive power output (MVAr)
Qgmax = gen(:, QMAX);   % gen max reactive power output (MVAr)

% -------- generalize def. of gen power output capacity to all bus -------- %
Pmin = Cg * Pgmin;
Pmax = Cg * Pgmax;
Qmin = Cg * Qgmin;
Qmax = Cg * Qgmax;

Pd = bus(:, PD);    % active power demand (MW)
Qd = bus(:, QD);    % reactive power demand (MVAr)

Vmin = bus(:, VMIN);    % min voltage magnitude (p.u.)
Vmax = bus(:, VMAX);    % max voltage magnitude (p.u.)

Smax = branch(:, RATE_A);   % branch flow limit (MVA)

% ----------------- unit conversion ---------------- %
C2 = C2*(baseMVA^2);      % $/h
C1 = C1*baseMVA;          % $/h
Pmin = Pmin/baseMVA;      % p.u.
Pmax = Pmax/baseMVA;      % p.u.
Qmin = Qmin/baseMVA;      % p.u.
Qmax = Qmax/baseMVA;      % p.u.
Pd = Pd/baseMVA;          % p.u.
Qd = Qd/baseMVA;          % p.u.
Smax = Smax/baseMVA;      % p.u.

%% build system admittance matrix
[Yk, Yk_, Mk, Ylineft, Ylinetf, Ylineft_, Ylinetf_] = makesdpmat('pglib_opf_case3_lmbd');
emat = speye(nb);
e = @(k) emat(:,k);
Nk = @(k) blkdiag(zeros(nb,nb), e(k)*e(k).');

%% solve the SDP relaxation of the OPF
alpha = sdpvar(ng, 1, 'full');          % auxiliary variable for quadratic cost function ($/h)
W = sdpvar(2*nb, 2*nb, 'symmetric');    % 2nb * 2nb real symmetric matrix, i.e. 
                                        % [real(V) imag(V)]' * [real(V) imag(V)] (p.u.)

obj = sum(alpha);   % $/h
const = [];

% generator active power limit (p.u.)
for k = 1:nb
    const = [const, Pmin(k) - Pd(k) <= trace(Yk(k)*W) <= Pmax(k) - Pd(k)];
end

% generator reactive power limit (p.u.)
for k = 1:nb
    const = [const, Qmin(k) - Qd(k) <= trace(Yk_(k)*W) <= Qmax(k) - Qd(k)];
end

% voltage mangnitude limit (p.u.)
for k = 1:nb
    const = [const, Vmin(k)^2 <= trace(Mk(k)*W) <= Vmax(k)^2];
end

% apparent power flow limit measured at the "from" bus (p.u.)
for l = 1:nl
    const = [const, [   -Smax(l)^2             trace(Ylineft(l)*W) trace(Ylineft_(l)*W);
                     trace(Ylineft(l)*W)                -1                  0;
                     trace(Ylineft_(l)*W)                0                 -1] <= 0];
end

% apparent power flow measured at the "to" bus (p.u.)
for t = 1:nl
    const = [const, [   -Smax(t)^2             trace(Ylinetf(t)*W) trace(Ylinetf_(t)*W);
                     trace(Ylinetf(t)*W)                -1                  0;
                     trace(Ylinetf_(t)*W)                0                 -1] <= 0];
end

% gen cost constraint f_k(Pg_k) <= alpha_k (p.u.)
for g = 1:ng
    const = [const, [C1(g)*(trace(Yk(g)*W) + Pd(g)) + C0(g) - alpha(g)   sqrt(C2(g))*(trace(Yk(g)*W) + Pd(g));
                            sqrt(C2(g))*(trace(Yk(g)*W) + Pd(g))                      -1                    ] <= 0];
end

% positive semi-definite constraint
const = [const, W >= 0];

% referenc bus
const = [const, trace(Nk(1)*W) == 0];

opts = sdpsettings('solver', 'mosek', 'verbose', 2);
diags = optimize(const, obj, opts);

% ------------------------ post-processing ------------------------ %
% power injection
for k = 1:nb
    Pinj(k,1) = trace(Yk(k)*value(W));  % p.u.
    Qinj(k,1) = trace(Yk_(k)*value(W)); % p.u.
end

% apparent power flow
for l = 1:nl
    Sf(l,1) = abs(trace(Ylineft(l)*value(W)) + 1i*trace(Ylineft_(l)*value(W))); % p.u.
    St(l,1) = abs(trace(Ylinetf(l)*value(W)) + 1i*trace(Ylinetf_(l)*value(W))); % p.u.
end

[V, D] = eig(value(W));
V = V * diag(1./sqrt(sum(V.^2))');  % normalize the column of V
neigv = nnz(diag(D) >= 1e-2);       % number of nonzereo eigvalue of W
if neigv > 2
    disp('physically menaingless solution');
else
    [~, idx] = max(diag(D));
    rst.Vopt = sqrt(D(idx,idx)) * (V(1:nb,idx) + 1i*V(nb+1:end,idx)); % p.u.
    rst.Vm = abs(rst.Vopt); % p.u.
    rst.Va = rad2deg(angle(rst.Vopt));  % degree
    rst.Pg = Cg'*(Pinj+Pd)*baseMVA;   % MW
    rst.Qg = Cg'*(Qinj+Qd)*baseMVA;   % MVAr
    rst.Sf = Sf*baseMVA;  % MVA
    rst.St = St*baseMVA;  % MVA
end

%% solve the dual problem of the opf
lambdaLB = sdpvar(nb, 1, 'full');   % LM assoc. with the lower bound of gen active power ($/h)
lambdaUB = sdpvar(nb, 1, 'full');   % LM assoc. with the upper bound of gen active power ($/h)
gammaLB = sdpvar(nb, 1, 'full');    % LM assoc. with the lower bound of gen reactive power ($/h)
gammaUB = sdpvar(nb, 1, 'full');    % LM assoc. with the upper bound of gen reactive power ($/h)
muLB = sdpvar(nb, 1, 'full');       % LM assoc. with the lower bound of voltage magnitude ($/h)
muUB = sdpvar(nb, 1, 'full');       % LM assoc. with the upper bound of voltage magnitude ($/h)
Hft1 = sdpvar(3, 3, 'symmetric');   % GLM assoc. with the apparent power flow limit ($/h)
Hft2 = sdpvar(3, 3, 'symmetric');   % GLM assoc. with the apparent power flow limit ($/h)
Hft3 = sdpvar(3, 3, 'symmetric');   % GLM assoc. with the apparent power flow limit ($/h)
Htf1 = sdpvar(3, 3, 'symmetric');   % GLM assoc. with the apparent power flow limit ($/h)
Htf2 = sdpvar(3, 3, 'symmetric');   % GLM assoc. with the apparent power flow limit ($/h)
Htf3 = sdpvar(3, 3, 'symmetric');   % GLM assoc. with the apparent power flow limit ($/h)
R1 = sdpvar(2, 2, 'symmetric');     % GLM assoc. with the gen cost constraint ($/h)
R2 = sdpvar(2, 2, 'symmetric');     % GLM assoc. with the gen cost constraint ($/h)
R3 = sdpvar(2, 2, 'symmetric');     % GLM assoc. with the gen cost constraint ($/h)

% ---------------- aggregate Lagrange multipliers and functions ---------------- %
lambda = lambdaUB - lambdaLB + Cg*(C1 + 2*sqrt(C2).*[R1(1,2);R2(1,2);R3(1,2)]);
gamma = gammaUB - gammaLB;
mu = muUB - muLB;

h = sum(lambdaLB.*Pmin - lambdaUB.*Pmax + lambda.*Pd + gammaLB.*Qmin - gammaUB.*Qmax + ...
    gamma.*Qd + muLB.*(Vmin.^2) - muUB.*(Vmax.^2), 'all') + ...
    sum(C0 - [R1(2,2);R2(2,2);R3(2,2)], 'all') - ...
    sum([Smax.^2;Smax.^2].*[Hft1(1,1);Hft2(1,1);Hft3(1,1);Htf1(1,1);Htf2(1,1);Htf3(1,1)] + ...
    [Hft1(2,2);Hft2(2,2);Hft3(2,2);Htf1(2,2);Htf2(2,2);Htf3(2,2)] + ...
    [Hft1(3,3);Hft2(3,3);Hft3(3,3);Htf1(3,3);Htf2(3,3);Htf3(3,3)], 'all');
A = 0;
for k = 1:nb
    A = A + lambda(k)*Yk(k) + gamma(k)*Yk_(k) + mu(k)*Mk(k);
end
A = A + 2*Hft1(1,2)*Ylineft(1) + 2*Hft1(1,3)*Ylineft_(1) + ...
        2*Hft2(1,2)*Ylineft(2) + 2*Hft2(1,3)*Ylineft_(2) + ...
        2*Hft3(1,2)*Ylineft(3) + 2*Hft3(1,3)*Ylineft_(3) + ...
        2*Htf1(1,2)*Ylinetf(1) + 2*Htf1(1,3)*Ylinetf_(1) + ...
        2*Htf2(1,2)*Ylinetf(2) + 2*Htf2(1,3)*Ylinetf_(2) + ...
        2*Htf3(1,2)*Ylinetf(3) + 2*Htf3(1,3)*Ylinetf_(3);

obj = -h;   % maximize dual function ($/h)
const = [];
const = [const, A >= 0];
const = [const, Hft1 >= 0, Hft2 >= 0, Hft3 >= 0, Htf1 >= 0, Htf2 >= 0, Htf3 >= 0];
const = [const, R1 >= 0, R1(1,1) == 1, R2 >= 0, R2(1,1) == 1, R3 >= 0, R3(1,1) == 1];
const = [const, lambdaLB >= 0, lambdaUB >= 0, gammaLB >= 0, gammaUB >= 0, muLB >= 0, muUB >= 0];

opts = sdpsettings('solver', 'mosek', 'verbose', 2);
diags = optimize(const, obj, opts);

% ------------------------ post-processing ------------------------ %
rst.lambda = value(lambda)/baseMVA;   % $/MWh
rst.gamma = value(gamma)/baseMVA;     % $/MVArh
rst.cost = -value(obj);               % $/h
rst.dimN_A = nnz(eig(value(A)) <1e-1);

disp('done!');