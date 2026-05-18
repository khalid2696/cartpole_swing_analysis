% =========================================================================
%
% Ellipsoidal tube reachability analysis for the cart-pole swing-up.
% Given a nominal trajectory (x_nom, u_nom) from direct collocation and
% a TVLQR struct, this script:
%   1. Interpolates S(t), K(t), x_nom(t), u_nom(t) onto a fine grid
%   2. Linearises the dynamics numerically at each grid point
%   3. Computes alpha(t) and beta(t) for the tube ODE
%   4. Integrates rho(t) forward
%   5. Checks ellipsoidal containment at t_f
%
% INPUTS (expected in workspace before running):
%   tvlqr       - struct with fields: t, S, K, Q, R
%                   t   : [1 x N]     time knot points
%                   S   : [4 x 4 x N] Riccati matrix at knot points
%                   K   : [1 x 4 x N] feedback gain at knot points
%                   Q   : [4 x 4]     state cost matrix
%                   R   : [1 x 1]     input cost
%   x_nom       : [4 x N]  nominal state trajectory at knot points
%   u_nom       : [1 x N]  nominal control at knot points
%   params      : struct with fields M, m, L, g, F_max
%   eps0        : scalar  -- size of initial ellipsoid X0 (rho at t=0)
%   eps_f       : scalar  -- size of terminal ellipsoid Xf
%   P_f         : [4 x 4] -- shape matrix of Xf (e.g. eye(4) for sphere)
%
% =========================================================================
clc; clearvars; close all
addpath("./lib/");
%% ── Initialization: nominal trajectory and TVLQR feedback gains

% Load the nominal trajectory and feedback control, and 
load('./precomputedData/swing_down/nominal_trajectory_and_input.mat');
load('./precomputedData/swing_down/TVLQR_gains_and_cost_matrices.mat');

%% ── Inputs ───────────────────────────────────────────────────

eps0  = 1e-4;        % initial ellipsoid size -- start small, grow until it fails
eps_f = 1/192.2777;         % terminal ellipsoid size -- your choice based on X_f spec
P_f   = eye(4);      % shape of X_f -- replace with your actual terminal set matrix

%% ── 0. Fine grid setup ───────────────────────────────────────────────────

N_fine   = 1000;                          % number of fine grid points
t_knots  = tvlqr.t;                       % original knot points  [1 x N]
t0       = t_knots(1);
tf       = t_knots(end);
t_fine   = linspace(t0, tf, N_fine);      % uniform fine grid      [1 x N_fine]

fprintf('Time horizon: [%.3f, %.3f] s,  N_fine = %d\n', t0, tf, N_fine);

%% ── 1. Interpolate S(t), K(t), x_nom(t), u_nom(t) ───────────────────────
% S is [4x4xN], K is [1x4xN] -- reshape to 2-D for interp1 then reshape back

nS = size(tvlqr.S, 1);   % = 4
nK = size(tvlqr.K, 1);   % = 1 (scalar input)

% Flatten trailing dimension for interpolation: [N x (4*4)] etc.
S_flat = reshape(tvlqr.S, nS*nS, [])';     % [N x 16]
K_flat = reshape(tvlqr.K, nK*nS, [])';     % [N x 4]
xn_T   = x_nom';                            % [N x 4]
un_T   = u_nom';                            % [N x 1]

% Interpolate onto fine grid (pchip preserves positivity better than spline)
S_fine_flat = interp1(t_knots, S_flat, t_fine, 'pchip');   % [N_fine x 16]
K_fine_flat = interp1(t_knots, K_flat, t_fine, 'pchip');   % [N_fine x 4]
xnom_fine   = interp1(t_knots, xn_T,   t_fine, 'pchip')';  % [4 x N_fine]
unom_fine   = interp1(t_knots, un_T,   t_fine, 'pchip')';  % [1 x N_fine]

% Reshape back to [4 x 4 x N_fine] and [1 x 4 x N_fine]
S_fine = reshape(S_fine_flat', nS, nS, N_fine);  % [4 x 4 x N_fine]
K_fine = reshape(K_fine_flat', nK, nS, N_fine);  % [1 x 4 x N_fine]

fprintf('Interpolation done.\n');

%% ── 2. Numerical linearisation A(t), B(t) ───────────────────────────────
% Compute A(t) = df/dx and B(t) = df/du at each fine grid point via
% central finite differences.

eps_fd = 1e-6;   % finite difference step

A_fine = zeros(4, 4, N_fine);  % [4 x 4 x N_fine]
B_fine = zeros(4, 1, N_fine);  % [4 x 1 x N_fine]

for k = 1:N_fine
    xk = xnom_fine(:, k);
    uk = unom_fine(k);

    % ── df/dx  (central differences, 4 columns) ──
    for i = 1:4
        ei          = zeros(4,1);
        ei(i)       = eps_fd;
        A_fine(:,i,k) = (cartPoleDynamics(xk+ei, uk, params) ...
                       - cartPoleDynamics(xk-ei, uk, params)) / (2*eps_fd);
    end

    % ── df/du  (central difference, scalar input) ──
    B_fine(:,1,k) = (cartPoleDynamics(xk, uk+eps_fd, params) ...
                   - cartPoleDynamics(xk, uk-eps_fd, params)) / (2*eps_fd);
end

fprintf('Linearisation done.\n');

%% ── 3. Closed-loop matrix A_cl(t) = A(t) - B(t)*K(t) ───────────────────

Acl_fine = zeros(4, 4, N_fine);
for k = 1:N_fine
    Acl_fine(:,:,k) = A_fine(:,:,k) ...
                    - B_fine(:,:,k) * K_fine(:,:,k);
end

%% ── 4. Compute alpha(t) and beta(t) ─────────────────────────────────────
%
%  alpha(t) = lambda_min(Q) / lambda_max(S(t))
%       -- LTV Lyapunov decay rate
%
%  beta(t)  = conservative Lipschitz bound on the nonlinear residual.
%             We use:
%               beta(t) = 2 * ||S(t)||_2 * L2(t) / lambda_min(S(t))
%             where L2(t) is an upper bound on the second derivative of f
%             w.r.t. x along the nominal, estimated via finite differences
%             of A(t) w.r.t. x (i.e., the "Jacobian of the Jacobian").
%
%  Derivation sketch for beta:
%    |nabla_V' g(dx,t)| = |2 dx' S g|
%    <= 2 ||S||_2 ||dx|| ||g||
%    g = O(||dx||^2)  =>  ||g|| <= (L2/2) ||dx||^2  on E(rho)
%    On E(rho): ||dx||^2 <= rho / lambda_min(S)
%    So |nabla_V' g| <= 2 ||S|| (L2/2) (rho/lmin_S) ||dx||
%                    <= ||S|| L2 / lmin_S * (dx' S dx)  [using ||dx|| <= sqrt(rho/lmin_S)]
%    => beta(t) = ||S(t)||_2 * L2(t) / lambda_min(S(t))
%
%  L2(t) is estimated as max_i ||dA_i/dx||_F  where A_i is the i-th row
%  of A(t,x), evaluated via a second finite difference in x.

Q        = tvlqr.Q;
lmin_Q   = min(eig(Q));   % scalar, constant

alpha_fine = zeros(1, N_fine);
beta_fine  = zeros(1, N_fine);
lmin_S     = zeros(1, N_fine);

eps_fd2 = 1e-4;   % larger step for second-order FD (A changes slowly)

for k = 1:N_fine
    Sk  = S_fine(:,:,k);
    xk  = xnom_fine(:, k);
    uk  = unom_fine(k);

    ev_S      = eig(Sk);
    lmax_S_k  = max(ev_S);
    lmin_S_k  = min(ev_S);
    lmin_S(k) = lmin_S_k;

    % alpha(t)
    alpha_fine(k) = lmin_Q / lmax_S_k;

    % L2(t): estimate ||d^2 f / dx^2|| via second FD of f w.r.t. x
    % Use the Frobenius norm of the finite-difference Hessian row-by-row
    L2_k = 0;
    for i = 1:4
        ei = zeros(4,1); ei(i) = eps_fd2;
        % Second derivative of f w.r.t. x_i x_i (diagonal of Hessian)
        d2f = (cartPoleDynamics(xk+ei, uk, params) ...
             - 2*cartPoleDynamics(xk,    uk, params) ...
             + cartPoleDynamics(xk-ei,  uk, params)) / eps_fd2^2;
        L2_k = max(L2_k, norm(d2f));
    end

    % beta(t)
    normS_k       = norm(Sk, 2);   % spectral norm = lambda_max(S)
    beta_fine(k)  = normS_k * L2_k / lmin_S_k;
end

fprintf('alpha and beta computed.\n');
fprintf('  mean(alpha) = %.4f,  max(beta) = %.4f\n', ...
        mean(alpha_fine), max(beta_fine));

%% ── 5. Integrate rho(t) forward ─────────────────────────────────────────
%
%  rho_dot = (beta(t) - alpha(t)) * rho(t),    rho(0) = rho0
%
%  We use ode45 with the fine-grid (alpha-beta) as an interpolated input.

rho0 = eps0;   % initial ellipsoid size = size of X0

% Build interpolant for (beta - alpha) on fine grid
gamma_fine = beta_fine - alpha_fine;   % [1 x N_fine]

gamma_interp = @(t) interp1(t_fine, gamma_fine, t, 'pchip', 'extrap');
rho_ode      = @(t, rho) gamma_interp(t) .* rho;

[t_rho, rho_sol] = ode45(rho_ode, [t0, tf], rho0);

% Interpolate rho back onto fine grid for plotting / checking
rho_fine = interp1(t_rho, rho_sol, t_fine, 'pchip');

fprintf('rho integration done.\n');
fprintf('  rho(0)  = %.6f\n', rho_fine(1));
fprintf('  rho(tf) = %.6f\n', rho_fine(end));

%% ── 6. Containment check at t_f ─────────────────────────────────────────
%
%  E(tf; rho(tf)) = { dx : dx' S(tf) dx <= rho(tf) }
%  Xf             = { x  : (x-xf)' Pf^{-1} (x-xf) <= eps_f }
%                 = { dx : dx' Pf^{-1} dx <= eps_f }   (in deviation coords)
%
%  Containment: E(tf; rho(tf)) ⊆ Xf
%  iff  rho(tf) * Pf^{-1} ⪯ eps_f * S(tf)
%  iff  eps_f * S(tf) - rho(tf) * Pf^{-1} is positive semidefinite.
%
%  Equivalently (avoid inverting Pf):
%  iff  rho(tf) / eps_f * S(tf)^{-1} ⪯ Pf
%  iff  all eigenvalues of  (rho(tf)/eps_f) * S(tf)^{-1} - Pf  are <= 0.

Stf     = S_fine(:,:,end);
rho_tf  = rho_fine(end);

% Build the matrix to check: eps_f * S(tf) - rho(tf) * inv(P_f)
M_check = eps_f * Stf - rho_tf * inv(P_f);
ev_check = eig(M_check);
certified = all(ev_check >= -1e-8);   % allow small numerical tolerance

fprintf('\n─── Containment Check ───────────────────────────────────\n');
fprintf('  rho(tf)              = %.6f\n', rho_tf);
fprintf('  eps_f                = %.6f\n', eps_f);
fprintf('  min eig(M_check)     = %.6e\n', min(ev_check));
if certified
    fprintf('  RESULT: ✓  E(tf;rho(tf)) ⊆ Xf  -- reachability CERTIFIED\n');
else
    fprintf('  RESULT: ✗  Tube does not fit inside Xf at tf.\n');
    fprintf('             Consider: smaller eps0, longer tf, or SOS method.\n');
end
fprintf('─────────────────────────────────────────────────────────\n\n');

%% ── 7. Plots ─────────────────────────────────────────────────────────────

figure('Name','Riccati Tube Analysis','NumberTitle','off');

% ── 7a. rho(t) ──
subplot(2,2,1);
plot(t_fine, rho_fine, 'b-', 'LineWidth', 1.5); hold on;
yline(eps_f, 'r--', '\epsilon_f', 'LabelHorizontalAlignment','left');
xlabel('Time (s)'); ylabel('\rho(t)');
title('Ellipsoid Scale \rho(t)');
legend('\rho(t)','\epsilon_f','Location','northwest');
grid on;

% ── 7b. alpha(t) and beta(t) ──
subplot(2,2,2);
plot(t_fine, alpha_fine, 'g-',  'LineWidth', 1.5); hold on;
plot(t_fine, beta_fine,  'r-',  'LineWidth', 1.5);
plot(t_fine, gamma_fine, 'k--', 'LineWidth', 1.0);
xlabel('Time (s)');
legend('\alpha(t)','\beta(t)','\beta-\alpha','Location','best');
title('Decay (\alpha) and Residual Bound (\beta)');
yline(0, 'k:', 'HandleVisibility','off');
grid on;

% ── 7c. Tube width along each state dimension ──
% Tube half-width in state i: sqrt(rho(t) / (e_i' S(t) e_i))
subplot(2,2,3);
state_labels = {'x_c (m)', '\dot{x}_c (m/s)', '\theta (rad)', '\dot\theta (rad/s)'};
colors = lines(4);
for i = 1:4
    ei = zeros(4,1); ei(i) = 1;
    half_width = zeros(1, N_fine);
    for k = 1:N_fine
        Sk = S_fine(:,:,k);
        half_width(k) = sqrt(rho_fine(k) / (ei' * Sk * ei));
    end
    plot(t_fine, half_width, 'Color', colors(i,:), 'LineWidth', 1.5); hold on;
end
xlabel('Time (s)'); ylabel('Half-width');
title('Tube Half-Width per State Dimension');
legend(state_labels, 'Location','best');
grid on;

% ── 7d. lambda_min(S(t)) ──
subplot(2,2,4);
plot(t_fine, lmin_S, 'm-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('\lambda_{min}(S(t))');
title('Minimum Eigenvalue of S(t)');
grid on;

sgtitle('Ellipsoidal Tube Reachability Analysis', 'FontSize', 13);

%% ── Helper: cart-pole dynamics ──────────────────────────────────────────
% Place this in a separate file cartPoleDynamics.m or as a local function.

function f = cartPoleDynamics(x, u, p)
% cartPoleDynamics  Evaluate the cart-pole vector field f(x,u).
%
%   x = [x_cart; v_x; theta; omega]
%   u = horizontal force (scalar)
%   p = struct with fields M, m, L, g

    x_cart = x(1);
    v_x    = x(2);
    theta  = x(3);
    omega  = x(4);

    s = sin(theta);
    c = cos(theta);
    D = p.M + p.m * s^2;

    f = [
        v_x;
        (u + p.m*p.L*omega^2*s + p.m*p.g*s*c) / D;
        omega;
        (-u*c - p.m*p.L*omega^2*s*c - (p.M+p.m)*p.g*s) / (p.L * D)
    ];
end