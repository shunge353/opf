function [Yk, Yk_, Mk, Ylineft, Ylinetf, Ylineft_, Ylinetf_] = makesdpmat(mpc)
% Create functions that return the matrices used in the semidefinite
% programming relaxation of the optimal power flow problem.
%
%   Inputs:
%       mpc: MATPOWER case variable with internal indexing
%   Outputs:
%       Yk: Function to create the matrix for active power injection.
%           Yk(i) is the matrix corresponding to the active power injection
%           at bus i.
%       Yk_: Function to create the matrix for reactive power injection.
%           Yk_(i) is the matrix corresponding to the reactive power
%           injection at bus i.
%       Mk: Function to create the matrix for the square of voltage
%           magnnitude. Mk(i) is the matrix corresponding to the square of
%           the voltage magnitude at bus i.
%       Ylineft: Function to create the matrix for the active power flow
%           on the specified line, measured from the "from" bus to the "to"
%           bus. Ylineft(i) is the matrix corresponding to the active power
%           flow on the line mpc.branch(i,:).
%       Ylinetf: Function to create the matrix for the reactive power flow
%           on the specified line, measured from the "to" bus to the "from"
%           bus. Ylinetf(i) is the matrix corresponding to the active power
%           flow on the line mpc.branch(i,:).
%       Ylinft_: Function to create the matrix for the reactive power
%           flow on the specified line, measured from the :from" bus to the
%           "to" bus. Ylineft_(i) is the matrix corresponding to the
%           reactive power flow on the line mpc.branch(i,:).
%       Ylinetf_: Function to create the matrix for the reactive power
%           flow on the specified line, measured from the "to" bus to the
%           "from" bus. Ylinetf_(i) is the matrix corresponding to the
%           reactive power flow on the line mpc.branch(i,:).

%% setup
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

mpc = loadcase(mpc);
nbus = size(mpc.bus, 1);
Y = makeYbus(mpc);

emat = speye(nbus);
e = @(k) emat(:, k);    % kth standard basis vector
Yk_small = @(k) e(k)*e(k).'*Y;

%% matrices used in SDP relaxation of OPF problem
Yk = @(k) (1/2) * [real(Yk_small(k) + Yk_small(k).') imag(Yk_small(k).' - Yk_small(k));
    imag(Yk_small(k) - Yk_small(k).') real(Yk_small(k) + Yk_small(k).')];
% Yk = @(k) (1/2)*[real(Yk_small(k) + Yk_small(k).') imag(Yk_small(k).' - Yk_small(k));
%     imag(Yk_small(k) - Yk_small(k).') real(Yk_small(k) + Yk_small(k).')];

Yk_ = @(k) -(1/2) * [imag(Yk_small(k) + Yk_small(k).') real(Yk_small(k) - Yk_small(k).');
    real(Yk_small(k).' - Yk_small(k)) imag(Yk_small(k) + Yk_small(k).')];

Mk = @(k) blkdiag(e(k)*e(k).', e(k)*e(k).');

% For the line limit matrices, specify a line index corresponding to the
%   entries in mpc.branch
gl = @(lineidx) real( 1 / (mpc.branch(lineidx,BR_R) + 1i*mpc.branch(lineidx,BR_X)) ); % real part of line admittance
bl = @(lineidx) imag( 1 / (mpc.branch(lineidx,BR_R) + 1i*mpc.branch(lineidx,BR_X)) ); % imaginary part of line admittance
bsl = @(lineidx) mpc.branch(lineidx,BR_B); % line shunt susceptance

Ylineft_small = @(lineidx) ((1/2)*1i*bsl(lineidx) + gl(lineidx) + 1i*bl(lineidx)) * e(mpc.branch(lineidx,F_BUS)) * e(mpc.branch(lineidx,F_BUS)).' - ...
    (gl(lineidx) + 1i*bl(lineidx)) * e(mpc.branch(lineidx,F_BUS)) * e(mpc.branch(lineidx,T_BUS)).';
Ylinetf_small = @(lineidx) ((1/2)*1i*bsl(lineidx) + gl(lineidx) + 1i*bl(lineidx)) * e(mpc.branch(lineidx,T_BUS)) * e(mpc.branch(lineidx,T_BUS)).' - ...
    (gl(lineidx) + 1i*bl(lineidx)) * e(mpc.branch(lineidx,T_BUS)) * e(mpc.branch(lineidx,F_BUS)).';

Ylineft = @(lidx) (1/2) * [real(Ylineft_small(lidx) + Ylineft_small(lidx).') imag(Ylineft_small(lidx).' - Ylineft_small(lidx));
    imag(Ylineft_small(lidx) - Ylineft_small(lidx).') real(Ylineft_small(lidx) + Ylineft_small(lidx).')];

Ylinetf = @(lidx) (1/2) * [real(Ylinetf_small(lidx) + Ylinetf_small(lidx).') imag(Ylinetf_small(lidx).' - Ylinetf_small(lidx));
    imag(Ylinetf_small(lidx) - Ylinetf_small(lidx).') real(Ylinetf_small(lidx) + Ylinetf_small(lidx).')];

Ylineft_ = @(lidx) -(1/2) * [imag(Ylineft_small(lidx) + Ylineft_small(lidx).') real(Ylineft_small(lidx) - Ylineft_small(lidx).');
    real(Ylineft_small(lidx).' - Ylineft_small(lidx)) imag(Ylineft_small(lidx) + Ylineft_small(lidx).')];

Ylinetf_ = @(lidx) -(1/2) * [imag(Ylinetf_small(lidx) + Ylinetf_small(lidx).') real(Ylinetf_small(lidx) - Ylinetf_small(lidx).');
    real(Ylinetf_small(lidx).' - Ylinetf_small(lidx)) imag(Ylinetf_small(lidx) + Ylinetf_small(lidx).')];


end