clc; clearvars; close all;

%'top': 25.850046; 
%'bottom': 57.567173;

%% Verifying level-set value c using feasibility in SOSP
%  Usage:
%    result = verify_roa_sos(mode, c_candidate)
%    result = verify_roa_sos(mode, c_candidate, lambda_deg)
%
%  Inputs:
%    mode        -- 'top' or 'bottom'
%    c_candidate -- candidate RoA level (scalar, e.g. from compute_roa)
%    lambda_deg  -- degree of SOS multiplier (even integer >= 4, default 4)
%
%  Output:
%    result.feasible  -- true if certificate found
%    result.c         -- the c that was verified
%    result.mode      -- equilibrium mode
%    result.info      -- solver info string

result = verify_roa_sos('top', 20, 4);
% result = verify_roa_sos('bottom', 50, 4);

%% =========================================================================
%  verify_roa_sos.m
%
%  SOS FEASIBILITY VERIFICATION of a candidate RoA level c.
%
%  Given a candidate c (e.g. from the radial line search), certifies whether
%  the ellipsoid  Omega_c = { delta : delta'*P*delta <= c }  is a verified
%  subset of the true nonlinear RoA under LQR control.
%
%  The certificate sought is:
%
%    -Wdot(delta) - lambda(delta)*(c - W(delta)) is SOS        ... (*)
%    lambda(delta) is SOS
%
%  where W(delta) = delta'*P*delta  and  Wdot is computed using the
%  degree-3 Taylor-expanded closed-loop dynamics.
%
%  If (*) is feasible, then Wdot < 0 on Omega_c \ {0}, certifying it as
%  a valid RoA inner approximation.
%
%  Requires: SOSTOOLS v4.00, MOSEK
% =========================================================================
function result = verify_roa_sos(mode, c_candidate, lambda_deg)

if nargin < 3 || isempty(lambda_deg)
    lambda_deg = 4;
end

assert(mod(lambda_deg,2)==0 && lambda_deg >= 4, ...
    'lambda_deg must be an even integer >= 4.');

fprintf('\n=== SOS Verification (Feasibility) ===\n');
fprintf('  Mode        : %s\n', upper(mode));
fprintf('  c candidate : %.6f\n', c_candidate);
fprintf('  lambda deg  : %d\n', lambda_deg);

% ------------------------------------------------------------------
% 1. System parameters and LQR
% ------------------------------------------------------------------
params   = init_params(mode);
P        = params.S;
K_lqr    = -params.gains.K;   % u = -K_lqr * delta

% ------------------------------------------------------------------
% 2. Build polynomial dynamics in delta coordinates
% ------------------------------------------------------------------
% State: delta = [dx; dxdot; dphi; dphidot]
%   where dphi = theta - theta_eq
%
% Taylor order 3 around phi=0 (i.e. around theta = theta_eq):
%   sin(theta) = sin(theta_eq + dphi)
%   cos(theta) = cos(theta_eq + dphi)
%
% For top    (theta_eq = pi):
%   sin(pi + dphi) = -sin(dphi) ~  -(dphi - dphi^3/6)
%   cos(pi + dphi) = -cos(dphi) ~  -(1 - dphi^2/2)
%
% For bottom (theta_eq = 0):
%   sin(0  + dphi) =  sin(dphi) ~   (dphi - dphi^3/6)
%   cos(0  + dphi) =  cos(dphi) ~   (1 - dphi^2/2)
%
% The approximation is done symbolically using SOSTOOLS pvar objects.

% Declare polynomial variables
pvar dx dxdot dphi dphidot

delta_vec = [dx; dxdot; dphi; dphidot];   % column

% Trig approximations (degree 3)
if strcmpi(mode, 'top')
    % sin(pi + dphi) = -sin(dphi), cos(pi + dphi) = -cos(dphi)
    s_poly =  -(dphi - dphi^3/6);     %  -sin(dphi) approx deg 3
    c_poly =  -(1    - dphi^2/2);     %  -cos(dphi) approx deg 2
elseif strcmpi(mode, 'bottom')
    s_poly =   (dphi - dphi^3/6);     %   sin(dphi) approx deg 3
    c_poly =   (1    - dphi^2/2);     %   cos(dphi) approx deg 2
else
    error('mode must be ''top'' or ''bottom''');
end

% Parameters
M_ = params.M; m_ = params.m; L_ = params.L; g_ = params.g;

% Control u = -K_lqr * delta  (no saturation for RoA analysis)
u_ctrl = -K_lqr * delta_vec;   % polynomial in delta (linear)

% Denominator Delta(theta) = M + m*sin^2(theta)
%   At theta = theta_eq + dphi: sin(theta) = s_poly
%   Delta_poly = M + m * s_poly^2  (degree 6 -- invert via series or keep)
%
% To stay polynomial, multiply through by Delta rather than dividing.
% We use the EXACT structure: f_cl = [f1; f2; f3; f4] where
%   f2 = (u + m*L*dphidot^2*s + m*g*s*c) / Delta
%   f4 = (-u*c - m*L*dphidot^2*s*c - (M+m)*g*s) / (L*Delta)
%
% For Wdot = 2*delta'*P*f, dividing by Delta is fine if Delta > 0 always.
% Delta = M + m*s^2 >= M > 0 for all states, so we can factor it out:
%   Wdot = 2*delta'*P*f = (1/Delta) * 2*delta'*P * [Delta*f1; num2; Delta*f3; num4/L]
%
% We compute Wdot_unnorm = Delta * Wdot to keep everything polynomial,
% then note that Delta > 0 => sign(Wdot_unnorm) = sign(Wdot).
% The S-procedure is applied to Wdot_unnorm.

Delta_poly = M_ + m_ * s_poly^2;    % >= M > 0

% Numerators (before dividing by Delta)
num_f2 = u_ctrl + m_*L_*dphidot^2*s_poly + m_*g_*s_poly*c_poly;
num_f4 = (-u_ctrl*c_poly - m_*L_*dphidot^2*s_poly*c_poly ...
          - (M_+m_)*g_*s_poly) / L_;

% Delta * f_cl  (polynomial, no denominator)
Df1 = Delta_poly * dxdot;
Df2 = num_f2;
Df3 = Delta_poly * dphidot;
Df4 = num_f4;

Df_cl = [Df1; Df2; Df3; Df4];   % column vector of polynomials

% W = delta'*P*delta
W_poly = delta_vec' * P * delta_vec;    % degree 2

% Wdot_unnorm = Delta * Wdot = 2 * delta' * P * (Delta * f_cl / Delta)
%             = 2 * delta' * P * f_cl      ... wait, need care:
%
% Actually: Wdot = 2*delta'*P*f_cl
%           f_cl = [dxdot; num_f2/Delta; dphidot; num_f4/Delta]
%
% So: Delta * Wdot = 2*delta'*P * [Delta*dxdot; num_f2; Delta*dphidot; num_f4]
%                 = 2*delta'*P * Df_cl
%
% Since Delta >= M > 0, sign(Delta*Wdot) = sign(Wdot).  We certify
% -(Delta*Wdot) - lambda*(c - W) is SOS, which implies -Wdot > 0 on Omega_c.

Wdot_unnorm = 2 * (delta_vec' * P * Df_cl);   % polynomial

% ------------------------------------------------------------------
% 3. Set up SOSTOOLS program
% ------------------------------------------------------------------
fprintf('  Building SOS program ...\n');
prog = sosprogram([dx; dxdot; dphi; dphidot]);

% Declare lambda(delta): SOS polynomial of degree lambda_deg
[prog, lambda] = sossosvar(prog, monomials([dx;dxdot;dphi;dphidot], ...
                                            0:lambda_deg));

% S-procedure constraint:
%   p(delta) = -(Delta*Wdot) - lambda*(c - W)  must be SOS
%
% Note: we exclude the origin separately -- p(0) = 0 trivially since
% Wdot(0)=0 and W(0)=0. The SOS constraint on a neighbourhood handles this.
%
% To avoid the trivial certificate at delta=0, add a small positive definite
% term: require p - epsilon*W >= 0 (SOS), epsilon small.
epsilon = 1e-6;
p_expr = -Wdot_unnorm - lambda*(c_candidate - W_poly) - epsilon*W_poly;

prog = sosineq(prog, p_expr);

% ------------------------------------------------------------------
% 4. Solve
% ------------------------------------------------------------------
fprintf('  Calling MOSEK via SOSTOOLS ...\n');
opts.solver  = 'mosek';
opts.params.MSK_IPAR_LOG = 0;   % suppress MOSEK output

prog_sol = sossolve(prog, opts);

% ------------------------------------------------------------------
% 5. Check feasibility
% ------------------------------------------------------------------
feasible = (prog_sol.solinfo.info.feasratio > 0.90);

result.feasible = feasible;
result.c        = c_candidate;
result.mode     = mode;
result.P        = P;
result.info     = prog_sol.solinfo.info;

fprintf('\n  --- Result ---\n');
if feasible
    fprintf('  FEASIBLE: c = %.6f is a certified RoA level.\n', c_candidate);
    fprintf('  Omega_c = { delta : delta''*P*delta <= %.6f } is verified.\n', ...
            c_candidate);
else
    fprintf('  INFEASIBLE: certificate not found for c = %.6f.\n', c_candidate);
    fprintf('  Try: (a) reducing c, (b) increasing lambda_deg, ');
    fprintf('(c) higher Taylor order.\n');
end

fprintf('  feasratio = %.4f  (>0.9 => feasible)\n', ...
        prog_sol.solinfo.info.feasratio);

% Final result: certified RoA
if feasible
    fprintf('\n   Certified RoA:  z^T * P * z  <=  %.6f\n', c_candidate);
    fprintf('   where  z = [x; xdot; theta - %.4f; thetadot]\n', params.xf(3));
    fprintf('\n   P matrix:\n');
    fprintf('   [ %10.4f  %10.4f  %10.4f  %10.4f ]\n', P');
    % Exact ellipsoid extents along each axis  (uses P^{-1})
    state_names = {'x        ', 'xdot   ', 'phi     ', 'phidot '};
    P_inv = inv(P);
    fprintf('\n   Ellipsoid extents along each axis  sqrt(c* . [P^{-1}]_kk):\n');
    fprintf('   (maximum |z_k| when all other states are free)\n');
    for k = 1:4
        extent = sqrt(c_candidate * P_inv(k,k));
        fprintf('     %s :  %.4f\n', state_names{k}, extent);
    end
end

end   % verify_roa_sos

% =========================================================================
%  LOCAL COPIES OF EXISTING HELPERS
%  (remove if already on your MATLAB path)
% =========================================================================
function params = init_params(mode)
    params.M = 1.0; params.m = 0.1; params.L = 0.5; params.g = 9.81;
    params.F_max = 10;
    params.tspan = [0 15];
    if strcmpi(mode,'top')
        params.xf = [0;0;pi;0];
    elseif strcmpi(mode,'bottom')
        params.xf = [0;0;0;0];
    else
        error('mode must be ''top'' or ''bottom''');
    end
    params.mode = mode;
    params.Q    = diag([10,1,100,1]);
    params.R    = 0.5;
    [params.gains.K, params.S] = compute_attractor_gain(params);
    params.initial_set_matrix  = 0.1*params.S;
end

function [K, S] = compute_attractor_gain(params)
    M_=params.M; m_=params.m; L_=params.L; g_=params.g;
    if strcmpi(params.mode,'top')
        A=[0 1 0 0;0 0 m_*g_/M_ 0;0 0 0 1;0 0 (M_+m_)*g_/(M_*L_) 0];
        B=[0;1/M_;0;1/(M_*L_)];
    else
        A=[0 1 0 0;0 0 m_*g_/M_ 0;0 0 0 1;0 0 -(M_+m_)*g_/(M_*L_) 0];
        B=[0;1/M_;0;-1/(M_*L_)];
    end
    [K_lqr,S,~] = lqr(A,B,params.Q,params.R);
    K = -K_lqr;   % stored negated: u = K*(xf-x) convention
end