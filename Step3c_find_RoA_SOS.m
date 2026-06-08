clc; clearvars; close all;

%'top': 25.850046; 
%'bottom': 57.567173;

%% Finding level-set value c using bisection search in SOSP

%  Usage:
%    result = find_roa_sos(mode, c_max)
%    result = find_roa_sos(mode, c_max, lambda_deg, tol, max_iter)
%
%  Inputs:
%    mode       -- 'top' or 'bottom'
%    c_max      -- upper bound on search (use c* from radial line search)
%    lambda_deg -- degree of SOS multiplier (even >= 4, default 4)
%    tol        -- bisection tolerance on c (default 1e-1)
%    max_iter   -- max bisection iterations (default 20)
%
%  Output:
%    result.c_sos      -- largest certified c found
%    result.c_max      -- upper bound supplied
%    result.feasible   -- true if at least one c > 0 was certified
%    result.history    -- table of bisection iterates [c_lo c_hi c_mid feas]
%    result.P          -- Riccati matrix
%    result.mode       -- equilibrium mode

% result = find_roa_sos('top', 26);
result = find_roa_sos('bottom', 58);

%% =========================================================================
%  find_roa_sos.m
%
%  SOS BISECTION to find the largest certifiable RoA level.
%
%  Searches for c_sos = max { c : SOS certificate exists for Omega_c }
%  by bisecting on c in [c_min, c_max], where c_max is supplied by the
%  caller (e.g. the radial line-search estimate c* from compute_roa).
%
%  At each bisection step, calls the same S-procedure feasibility check
%  as verify_roa_sos.  The bisection terminates when the interval width
%  falls below tol.
%
%  Certificate condition (S-procedure):
%    -(Delta*Wdot(delta)) - lambda(delta)*(c - W(delta)) is SOS
%    lambda(delta) is SOS
%  where Delta = M + m*sin^2(theta) > 0, absorbed to keep polynomials.
%
%  Requires: SOSTOOLS v4.00, MOSEK
% =========================================================================
function result = find_roa_sos(mode, c_max, lambda_deg, tol, max_iter)

% ---- defaults -----------------------------------------------------------
if nargin < 3 || isempty(lambda_deg), lambda_deg = 4;    end
if nargin < 4 || isempty(tol),        tol        = 1e-1; end
if nargin < 5 || isempty(max_iter),   max_iter   = 20;   end

assert(mod(lambda_deg,2)==0 && lambda_deg >= 4, ...
    'lambda_deg must be an even integer >= 4.');
assert(c_max > 0, 'c_max must be positive.');

fprintf('\n=== SOS Bisection (Optimization) ===\n');
fprintf('  Mode       : %s\n', upper(mode));
fprintf('  c_max      : %.6f  (from line search)\n', c_max);
fprintf('  lambda_deg : %d\n', lambda_deg);
fprintf('  tol        : %.2e\n', tol);
fprintf('  max_iter   : %d\n',   max_iter);

% ---- system / LQR -------------------------------------------------------
params  = init_params(mode);
P       = params.S;
K_lqr   = -params.gains.K;

% ---- build polynomial dynamics once (shared across all bisection steps) -
[Wdot_unnorm, W_poly, delta_vec] = build_polynomial_dynamics( ...
    mode, params, K_lqr);

% ---- bisection ----------------------------------------------------------
c_lo  = 0;
c_hi  = c_max;
c_sos = 0;          % best certified level so far

history = zeros(max_iter, 4);   % [c_lo, c_hi, c_mid, feasible]

fprintf('\n  Iter   c_lo        c_mid       c_hi        Feasible\n');
fprintf('  %-6s %-11s %-11s %-11s %s\n', ...
        '----','----------','----------','----------','--------');

for iter = 1:max_iter

    c_mid = (c_lo + c_hi) / 2;

    feas = sos_feasibility_check(c_mid, Wdot_unnorm, W_poly, ...
                                  delta_vec, lambda_deg);

    history(iter,:) = [c_lo, c_hi, c_mid, feas];

    fprintf('  %-6d %-11.5f %-11.5f %-11.5f %s\n', ...
            iter, c_lo, c_mid, c_hi, bool2str(feas));

    if feas
        c_sos = c_mid;   % best certified so far
        c_lo  = c_mid;   % search upper half
    else
        c_hi  = c_mid;   % search lower half
    end

    if (c_hi - c_lo) < tol
        fprintf('  Converged: interval width %.2e < tol (%.2e)\n', ...
                c_hi - c_lo, tol);
        break
    end
end

% ---- results ------------------------------------------------------------
feasible = (c_sos > 0);

result.c_sos    = c_sos;
result.c_max    = c_max;
result.feasible = feasible;
result.history  = history(1:iter, :);
result.P        = P;
result.mode     = mode;

fprintf('\n  --- Result ---\n');
if feasible
    fprintf('  Largest certified c_sos = %.6f\n', c_sos);
    fprintf('  (line-search c_max was  = %.6f)\n', c_max);
    ratio = c_sos / c_max;
    fprintf('  c_sos / c_max = %.4f\n', ratio);
    if ratio < 0.80
        fprintf('  NOTE: c_sos is substantially smaller than c_max.\n');
        fprintf('  Consider: (a) higher lambda_deg, (b) higher Taylor order.\n');
        fprintf('  The line-search c_max may be optimistic for this degree.\n');
    end
else
    fprintf('  No c > 0 could be certified.\n');
    fprintf('  Check: (a) lambda_deg >= 4, (b) Taylor approx valid,\n');
    fprintf('         (c) MOSEK license, (d) LQR gain is stabilising.\n');
end

% Final result: certified RoA
if feasible
    fprintf('\n   Certified RoA:  z^T * P * z  <=  %.6f\n', c_sos);
    fprintf('   where  z = [x; xdot; theta - %.4f; thetadot]\n', params.xf(3));
    fprintf('\n   P matrix:\n');
    fprintf('   [ %10.4f  %10.4f  %10.4f  %10.4f ]\n', P');
    % Exact ellipsoid extents along each axis  (uses P^{-1})
    state_names = {'x        ', 'xdot   ', 'phi     ', 'phidot '};
    P_inv = inv(P);
    fprintf('\n   Ellipsoid extents along each axis  sqrt(c* . [P^{-1}]_kk):\n');
    fprintf('   (maximum |z_k| when all other states are free)\n');
    for k = 1:4
        extent = sqrt(c_sos * P_inv(k,k));
        fprintf('     %s :  %.4f\n', state_names{k}, extent);
    end
end

end   % find_roa_sos

% =========================================================================
%  BUILD POLYNOMIAL DYNAMICS  (called once, reused across bisection steps)
% =========================================================================
function [Wdot_unnorm, W_poly, delta_vec] = build_polynomial_dynamics( ...
    mode, params, K_lqr)

pvar dx dxdot dphi dphidot
delta_vec = [dx; dxdot; dphi; dphidot];

% Degree-3 Taylor for sin/cos around theta_eq
if strcmpi(mode,'top')
    % theta = pi + dphi:  sin = -sin(dphi), cos = -cos(dphi)
    s_poly = -(dphi - dphi^3/6);
    c_poly = -(1    - dphi^2/2);
else
    % theta = 0  + dphi:  sin = sin(dphi),  cos = cos(dphi)
    s_poly =  (dphi - dphi^3/6);
    c_poly =  (1    - dphi^2/2);
end

M_=params.M; m_=params.m; L_=params.L; g_=params.g;
P = params.S;

u_ctrl     = -K_lqr * delta_vec;
Delta_poly = M_ + m_ * s_poly^2;

num_f2 = u_ctrl + m_*L_*dphidot^2*s_poly + m_*g_*s_poly*c_poly;
num_f4 = (-u_ctrl*c_poly - m_*L_*dphidot^2*s_poly*c_poly ...
          - (M_+m_)*g_*s_poly) / L_;

Df_cl = [Delta_poly*dxdot; num_f2; Delta_poly*dphidot; num_f4];

W_poly       = delta_vec' * P * delta_vec;
Wdot_unnorm  = 2 * (delta_vec' * P * Df_cl);

end

% =========================================================================
%  SOS FEASIBILITY CHECK  (called at each bisection step)
% =========================================================================
function feas = sos_feasibility_check(c, Wdot_unnorm, W_poly, ...
                                       delta_vec, lambda_deg)
epsilon = 1e-6;

prog = sosprogram(delta_vec);

[prog, lambda] = sossosvar(prog, monomials(delta_vec, 0:lambda_deg));

p_expr = -Wdot_unnorm - lambda*(c - W_poly) - epsilon*W_poly;
prog = sosineq(prog, p_expr);

opts.solver          = 'mosek';
opts.params.MSK_IPAR_LOG = 0;   % suppress solver output

try
    prog_sol = sossolve(prog, opts);
    feas = (prog_sol.solinfo.info.feasratio > 0.90);
catch
    feas = false;
end

end

% =========================================================================
%  UTILITY
% =========================================================================
function s = bool2str(b)
if b, s = 'YES'; else, s = 'NO'; end
end

% =========================================================================
%  LOCAL COPIES OF EXISTING HELPERS
%  (remove if already on your MATLAB path)
% =========================================================================
function params = init_params(mode)
    params.M=1.0; params.m=0.1; params.L=0.5; params.g=9.81;
    params.F_max=10; params.tspan=[0 15];
    if strcmpi(mode,'top')         
        params.xf=[0;0;pi;0];
    elseif strcmpi(mode,'bottom')  
        params.xf=[0;0;0;0];
    else 
        error('mode must be ''top'' or ''bottom'''); 
    end
    params.mode=mode;
    params.Q=diag([10,1,100,1]); params.R=0.5;
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
    [K_lqr,S,~]=lqr(A,B,params.Q,params.R);
    K=-K_lqr;
end