clc; clearvars; close all;
compute_roa();

% =========================================================================
%  compute_roa.m
%  Region of Attraction computation for cart-pole LQR controller
%  at the upright (top) and downward (bottom) equilibria.
%
%  Method (hybrid):
%    1. Radial line search  → accurate c* via fzero along random rays
%    2. 2-D slice grid      → visualisation of Wdot sign + W=c* contour
%
%  Uses existing helper functions:
%    init_params, compute_attractor_gain, cartpole_dynamics, LQR_attractor
%
%  Usage:
%    compute_roa          % runs both equilibria
%    compute_roa('top')   % upright only
%    compute_roa('bottom')% downward only
% =========================================================================
function compute_roa(mode_arg)

if nargin < 1
    modes = {'top', 'bottom'};
else
    modes = {mode_arg};
end

for mi = 1:numel(modes)
    mode = modes{mi};
    fprintf('\n=== RoA computation: %s equilibrium ===\n', upper(mode));

    % ------------------------------------------------------------------
    % 1. Initialise parameters and compute LQR gain / Riccati matrix
    % ------------------------------------------------------------------
    params = init_params(mode);
    P      = params.S;              % Riccati solution  (W = z'*P*z)
    K_lqr  = -params.gains.K;      % raw lqr K: u = -K_lqr*z
                                    % (params.gains.K stores -K_lqr for
                                    %  LQR_attractor compatibility)
    x_eq   = params.xf;            % equilibrium state [x; xdot; theta; thdot]

    fprintf('   LQR gain K_lqr = [%.4f  %.4f  %.4f  %.4f]\n', K_lqr);
    fprintf('   lambda(A-BK) = ');
    fprintf('%.4f  ', real(eig(compute_Acl(params))));
    fprintf('\n');

    % ------------------------------------------------------------------
    % 2.  Radial line search  →  c_star
    % ------------------------------------------------------------------
    fprintf('   Running radial line search ...\n');
    c_star = radial_line_search(P, K_lqr, x_eq, params);
    fprintf('   c* = %.6f\n', c_star);

    % ------------------------------------------------------------------
    % Print parametrised RoA  z'*P*z <= c*
    % ------------------------------------------------------------------
    fprintf('\n   Certified RoA:  z^T * P * z  <=  %.6f\n', c_star);
    fprintf('   where  z = [x; xdot; theta - %.4f; thetadot]\n', x_eq(3));
    fprintf('\n   P matrix:\n');
    fprintf('   [ %10.4f  %10.4f  %10.4f  %10.4f ]\n', P');
    % Exact ellipsoid extents along each axis  (uses P^{-1})
    state_names = {'x        ', 'xdot   ', 'phi     ', 'phidot '};
    P_inv = inv(P);
    fprintf('\n   Ellipsoid extents along each axis  sqrt(c* . [P^{-1}]_kk):\n');
    fprintf('   (maximum |z_k| when all other states are free)\n');
    for k = 1:4
        extent = sqrt(c_star * P_inv(k,k));
        fprintf('     %s :  %.4f\n', state_names{k}, extent);
    end
    
    % Axis-aligned extents  (all other states forced to zero — conservative)
    fprintf('\n [Conservative] Axis-aligned extents  sqrt(c*/P_kk):\n');
    fprintf('   (maximum |z_k| when all other states = 0)\n');
    for k = 1:4
        extent_diag = sqrt(c_star / P(k,k));
        fprintf('     %s :  %.4f\n', state_names{k}, extent_diag);
    end
    
    % ------------------------------------------------------------------
    % 3.  Physical validity check
    % ------------------------------------------------------------------
    %  Confirm Omega_{c*} does not extend past |phi| = pi/2.
    %  W(z) = z'Pz <= c*  and  z = [0,0,phi,0]
    %  => phi^2 * P(3,3) <= c*  => |phi| <= sqrt(c*/P(3,3))
    phi_max_roa = sqrt(c_star / P(3,3));
    fprintf('   Max |phi| in Omega_c* (other states=0): %.4f rad  (%.1f deg)\n', ...
            phi_max_roa, rad2deg(phi_max_roa));
    if strcmp(mode,'top') && phi_max_roa > pi/2
        warning('Omega_c* extends beyond pi/2 in phi — LQR validity questionable at boundary.');
    end

    % ------------------------------------------------------------------
    % 4.  2-D slice visualisation
    % ------------------------------------------------------------------
    fprintf('\n\n   Generating slice plots ...\n');
    plot_roa_slices(P, K_lqr, x_eq, c_star, params, mode);
end

end % main function

% =========================================================================
%  RADIAL LINE SEARCH
% =========================================================================
function c_star = radial_line_search(P, K, x_eq, params)
%  Shoot N_rays random unit directions in R^4.
%  Along each ray z(r) = r*d, find r* where Wdot first crosses zero.
%  c* = min W(z(r*)) over all rays that have a crossing.

N_rays   = 8000;    % number of random directions
r_max    = 8.0;     % maximum ray length (physical bound)
rng(42);            % reproducible random seed

% Random unit directions uniformly on S^3
D = randn(4, N_rays);
D = D ./ vecnorm(D, 2, 1);   % each column is a unit vector

c_candidates = inf(1, N_rays);

for i = 1:N_rays
    d = D(:, i);

    % Wdot at a scalar distance r along this ray
    wdot_fn = @(r) Wdot_scalar(r, d, P, K, x_eq, params);

    % Quick sign-check bracket:  must be negative near origin
    w0 = wdot_fn(1e-6);
    if w0 >= 0
        % Direction starts unstable (numerical noise near origin) — skip
        continue
    end

    % Find the far end of the bracket
    w_far = wdot_fn(r_max);
    if w_far < 0
        % Wdot stays negative all the way to r_max — not the limiting
        % direction, skip (contributes no finite c candidate)
        continue
    end

    % Binary-search bracket [0, r_max] has a sign change → use fzero
    try
        r_cross = fzero(wdot_fn, [1e-6, r_max], ...
                        optimset('TolX', 1e-8, 'Display', 'off'));
        z_cross  = r_cross * d;
        c_candidates(i) = z_cross' * P * z_cross;   % W at crossing
    catch
        % fzero failed for this ray — skip
    end
end

valid = c_candidates(isfinite(c_candidates));
if isempty(valid)
    error('No valid ray crossings found. Increase r_max or N_rays.');
end
c_star = min(valid);

fprintf('   Rays with crossings: %d / %d\n', numel(valid), N_rays);
fprintf('   Most restrictive c candidate: %.6f\n', c_star);
end

% =========================================================================
%  Wdot along a ray at scalar distance r
% =========================================================================
function wd = Wdot_scalar(r, d, P, K, x_eq, params)
% z = r*d,  q = x_eq + z,  u = -K*z  (raw lqr K, no saturation)
z    = r * d;
q    = x_eq + z;
u    = -K * z;            % standard stabilising law
dq   = cartpole_dynamics(0, q, u, params);
wd   = 2 * (z' * P * dq);
end

% =========================================================================
%  Wdot on a full state vector (for grid evaluation)
% =========================================================================
function wd = Wdot_state(z, P, K, x_eq, params)
q  = x_eq + z;
u  = -K * z;
dq = cartpole_dynamics(0, q, u, params);
wd = 2 * (z' * P * dq);
end

% =========================================================================
%  2-D SLICE PLOTS
% =========================================================================
function plot_roa_slices(P, K, x_eq, c_star, params, mode)
%  Three slices through the 4-D state space:
%    Slice A: (phi, phi_dot)   with  x=0, xdot=0
%    Slice B: (x,   xdot)     with  phi=0, phi_dot=0
%    Slice C: (x,   phi)      with  xdot=0, phi_dot=0
%
%  phi = theta - theta_eq

N_grid = 120;    % grid points per axis

% ---- axis ranges (physical) -------------------------------------------
if strcmp(mode, 'top')
    phi_rng  = linspace(-pi/2,  pi/2,  N_grid);
    om_rng   = linspace(-6,     6,     N_grid);
    x_rng    = linspace(-2,     2,     N_grid);
    xd_rng   = linspace(-4,     4,     N_grid);
    eq_label = 'Upright ($\theta=\pi$)';
else
    phi_rng  = linspace(-pi,    pi,    N_grid);
    om_rng   = linspace(-8,     8,     N_grid);
    x_rng    = linspace(-3,     3,     N_grid);
    xd_rng   = linspace(-5,     5,     N_grid);
    eq_label = 'Downward ($\theta=0$)';
end

% -----------------------------------------------------------------------
%  Build grids and evaluate Wdot + W for each slice
% -----------------------------------------------------------------------

% -- Slice A: (phi, omega) with x=0, xdot=0 ----------------------------
[PHI_A, OM_A] = meshgrid(phi_rng, om_rng);
W_A    = zeros(size(PHI_A));
Wd_A   = zeros(size(PHI_A));
for i = 1:N_grid
    for j = 1:N_grid
        z = [0; 0; PHI_A(i,j); OM_A(i,j)];
        W_A(i,j)  = z' * P * z;
        Wd_A(i,j) = Wdot_state(z, P, K, x_eq, params);
    end
end

% -- Slice B: (x, xdot) with phi=0, omega=0 ----------------------------
[X_B, XD_B] = meshgrid(x_rng, xd_rng);
W_B    = zeros(size(X_B));
Wd_B   = zeros(size(X_B));
for i = 1:N_grid
    for j = 1:N_grid
        z = [X_B(i,j); XD_B(i,j); 0; 0];
        W_B(i,j)  = z' * P * z;
        Wd_B(i,j) = Wdot_state(z, P, K, x_eq, params);
    end
end

% -- Slice C: (x, phi) with xdot=0, omega=0 ----------------------------
[X_C, PHI_C] = meshgrid(x_rng, phi_rng);
W_C    = zeros(size(X_C));
Wd_C   = zeros(size(X_C));
for i = 1:N_grid
    for j = 1:N_grid
        z = [X_C(i,j); 0; PHI_C(i,j); 0];
        W_C(i,j)  = z' * P * z;
        Wd_C(i,j) = Wdot_state(z, P, K, x_eq, params);
    end
end

% -----------------------------------------------------------------------
%  Figure
% -----------------------------------------------------------------------
fig = figure('Name', sprintf('RoA - %s equilibrium', mode), ...
             'Color', 'w', 'Position', [100 100 1300 420]);

% colour map: green = Wdot<0 (certified), red = Wdot>=0 (outside)
cmap_certified = [0.85 0.95 0.85;   % green tint (inside)
                  0.95 0.82 0.82];  % red  tint (outside)

% ── Slice A ─────────────────────────────────────────────────────────────
ax1 = subplot(1,3,1);
hold(ax1,'on');
% Shade by Wdot sign
pcolor(ax1, PHI_A, OM_A, double(Wd_A >= 0));
shading(ax1,'flat');
colormap(ax1, cmap_certified);
clim(ax1,[0 1]);
% W = c* contour
contour(ax1, PHI_A, OM_A, W_A, [c_star c_star], 'b-', 'LineWidth', 2.0);
% Wdot = 0 contour (black dashed)
contour(ax1, PHI_A, OM_A, Wd_A, [0 0], 'k--', 'LineWidth', 1.2);
plot(ax1, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
xlabel(ax1, '$\phi = \theta - \theta_{\rm eq}$ (rad)', ...
       'Interpreter','latex','FontSize',11);
ylabel(ax1, '$\dot\theta$ (rad s$^{-1}$)', ...
       'Interpreter','latex','FontSize',11);
title(ax1, sprintf('Slice A: $(\\phi,\\dot\\theta)$,  $x=\\dot x=0$\n%s', eq_label), ...
      'Interpreter','latex','FontSize',11);
grid(ax1,'on'); box(ax1,'on');
add_legend(ax1, c_star);

% ── Slice B ─────────────────────────────────────────────────────────────
ax2 = subplot(1,3,2);
hold(ax2,'on');
pcolor(ax2, X_B, XD_B, double(Wd_B >= 0));
shading(ax2,'flat');
colormap(ax2, cmap_certified);
clim(ax2,[0 1]);
contour(ax2, X_B, XD_B, W_B, [c_star c_star], 'b-', 'LineWidth', 2.0);
contour(ax2, X_B, XD_B, Wd_B, [0 0], 'k--', 'LineWidth', 1.2);
plot(ax2, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
xlabel(ax2, '$x$ (m)', 'Interpreter','latex','FontSize',11);
ylabel(ax2, '$\dot x$ (m s$^{-1}$)', 'Interpreter','latex','FontSize',11);
title(ax2, sprintf('Slice B: $(x,\\dot x)$,  $\\phi=\\dot\\theta=0$\n%s', eq_label), ...
      'Interpreter','latex','FontSize',11);
grid(ax2,'on'); box(ax2,'on');
add_legend(ax2, c_star);

% ── Slice C ─────────────────────────────────────────────────────────────
ax3 = subplot(1,3,3);
hold(ax3,'on');
pcolor(ax3, X_C, PHI_C, double(Wd_C >= 0));
shading(ax3,'flat');
colormap(ax3, cmap_certified);
clim(ax3,[0 1]);
contour(ax3, X_C, PHI_C, W_C, [c_star c_star], 'b-', 'LineWidth', 2.0);
contour(ax3, X_C, PHI_C, Wd_C, [0 0], 'k--', 'LineWidth', 1.2);
plot(ax3, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
xlabel(ax3, '$x$ (m)', 'Interpreter','latex','FontSize',11);
ylabel(ax3, '$\phi$ (rad)', 'Interpreter','latex','FontSize',11);
title(ax3, sprintf('Slice C: $(x,\\phi)$,  $\\dot x=\\dot\\theta=0$\n%s', eq_label), ...
      'Interpreter','latex','FontSize',11);
grid(ax3,'on'); box(ax3,'on');
add_legend(ax3, c_star);

% ── Super-title ──────────────────────────────────────────────────────────
sgtitle(fig, ...
    sprintf('LQR Region of Attraction - %s equilibrium   (c* = %.4f)', ...
            upper(mode), c_star), ...
    'Interpreter','latex','FontSize',13,'FontWeight','bold');

% ── Save ─────────────────────────────────────────────────────────────────
fname = sprintf('roa_%s', mode);
%saveas(fig, fullfile(fileparts(mfilename('fullpath')), [fname '.fig']));
print(fig, fullfile(fileparts(mfilename('fullpath')), [fname '.png']), ...
      '-dpng', '-r300');
fprintf('   Figure saved: %s.png\n', fname);
end

% =========================================================================
%  Legend helper
% =========================================================================
function add_legend(ax, c_star)
h1 = patch(ax, NaN, NaN, [0.85 0.95 0.85], 'EdgeColor','none');
h2 = patch(ax, NaN, NaN, [0.95 0.82 0.82], 'EdgeColor','none');
h3 = plot(ax, NaN, NaN, 'b-', 'LineWidth', 2.0);
h4 = plot(ax, NaN, NaN, 'k--','LineWidth', 1.2);
legend(ax, [h1 h2 h3 h4], ...
    {'$\dot W < 0$ (certified)', ...
     '$\dot W \geq 0$ (outside)', ...
     sprintf('$W = c^* = %.3f$', c_star), ...
     '$\dot W = 0$ boundary'}, ...
    'Interpreter','latex','FontSize',8,'Location','best');
end

% =========================================================================
%  Closed-loop A matrix  (for eigenvalue check)
% =========================================================================
function Acl = compute_Acl(params)
% params.gains.K = -K_lqr, so K_lqr = -params.gains.K
K_lqr = -params.gains.K;
M_ = params.M; m_ = params.m; L_ = params.L; g_ = params.g;
if strcmpi(params.mode, 'top')
    A = [0 1 0 0; 0 0 m_*g_/M_ 0; 0 0 0 1; ...
         0 0 (M_+m_)*g_/(M_*L_) 0];
    B = [0; 1/M_; 0; 1/(M_*L_)];
else
    A = [0 1 0 0; 0 0 m_*g_/M_ 0; 0 0 0 1; ...
         0 0 -(M_+m_)*g_/(M_*L_) 0];
    B = [0; 1/M_; 0; -1/(M_*L_)];
end
Acl = A - B*K_lqr;
end

% =========================================================================
%  EXISTING FUNCTIONS (kept here for self-contained script)
%  If these already exist in your codebase, remove them from here.
% =========================================================================

function u = LQR_attractor(x, params)
    K = params.gains.K;
    u = K * (params.xf - x);
    u = max(min(u, params.F_max), -params.F_max);
end

function params = init_params(mode)
    params.M = 1.0;
    params.m = 0.1;
    params.L = 0.5;
    params.g = 9.81;
    params.F_max = 10;
    params.tspan = [0 15];

    if strcmpi(mode, 'top')
        params.xf = [0; 0; pi; 0];
    elseif strcmpi(mode, 'bottom')
        params.xf = [0; 0; 0; 0];
    else
        error('mode must be ''top'' or ''bottom''');
    end

    params.mode = mode;
    params.Q    = diag([10, 1, 100, 1]);
    params.R    = 0.5;

    [params.gains.K, params.S] = compute_attractor_gain(params);
    params.initial_set_matrix  = 0.1 * params.S;
end

function dx = cartpole_dynamics(t, x, u, params)
    s = sin(x(3)); c = cos(x(3));
    v = x(2); omega = x(4);
    denom = params.M + params.m * s^2;
    f2 = (u + params.m*params.L*omega^2*s + params.m*params.g*s*c) / denom;
    f4 = (-u*c - params.m*params.L*omega^2*s*c - ...
          (params.M+params.m)*params.g*s) / (params.L * denom);
    dx = [v; f2; omega; f4];
end

function [K, S] = compute_attractor_gain(params)
    M_ = params.M; m_ = params.m; L_ = params.L; g_ = params.g;
    if strcmpi(params.mode, 'top')
        A = [0 1 0 0; 0 0 m_*g_/M_ 0; 0 0 0 1; ...
             0 0 (M_+m_)*g_/(M_*L_) 0];
        B = [0; 1/M_; 0; 1/(M_*L_)];
    elseif strcmpi(params.mode, 'bottom')
        A = [0 1 0 0; 0 0 m_*g_/M_ 0; 0 0 0 1; ...
             0 0 -(M_+m_)*g_/(M_*L_) 0];
        B = [0; 1/M_; 0; -1/(M_*L_)];
    else
        error('mode must be ''top'' or ''bottom''');
    end
    [K_lqr, S, ~] = lqr(A, B, params.Q, params.R);
    % lqr() stabilising law:  u = -K_lqr * z
    % LQR_attractor uses:     u = K_stored*(xf-x) = K_stored*(-z)
    % So store K_stored = -K_lqr to keep both consistent.
    K = -K_lqr;
end