function tvlqr = computeTVLQR(t_nom, x_nom, u_nom, Q, R, params)
    % COMPUTETVLQR  Compute the TVLQR struct for the cart-pole swing-up.
    %
    %   tvlqr = computeTVLQR(t_nom, x_nom, u_nom, Q, R, params)
    %
    % ── Inputs ───────────────────────────────────────────────────────────────
    %   t_nom   [1 x N]    time knot points (Drake-style, monotone increasing)
    %   x_nom   [4 x N]    nominal state trajectory:
    %                        row 1 : x_cart  (m)
    %                        row 2 : v_x     (m/s)
    %                        row 3 : theta   (rad)
    %                        row 4 : omega   (rad/s)
    %   u_nom   [1 x N]    nominal control (horizontal force, N)
    %   Q       [4 x 4]    positive-definite state cost matrix
    %   R       [1 x 1]    positive scalar input cost
    %   params  struct     cart-pole physical parameters:
    %                        .M      cart mass      (kg)
    %                        .m      pendulum mass  (kg)
    %                        .L      pendulum length (m)
    %                        .g      gravity         (m/s^2)
    %
    % ── Output ───────────────────────────────────────────────────────────────
    %   tvlqr   struct with fields:
    %     .t    [1 x N]       time knot points  (same as t_nom)
    %     .S    [4 x 4 x N]   Riccati matrix S(t)  at each knot point
    %     .K    [1 x 4 x N]   optimal feedback gain K(t) = R^{-1} B(t)' S(t)
    %     .Q    [4 x 4]       state cost (stored for downstream use)
    %     .R    [1 x 1]       input cost (stored for downstream use)
    %
    % ── Method ───────────────────────────────────────────────────────────────
    %   The time-varying Riccati differential equation (RDE):
    %
    %     -dS/dt = A(t)'S + S A(t) - S B(t) R^{-1} B(t)' S + Q
    %      S(tf) = Q                         (terminal condition)
    %
    %   is integrated BACKWARD from tf to t0.  We implement this by reversing
    %   the time axis: let tau = tf - t, so d/dtau = -d/dt, giving
    %
    %      dS/dtau = A(tf-tau)'S + S A(tf-tau) - S B(tf-tau) R^{-1} B(tf-tau)' S + Q
    %      S(tau=0) = Q
    %
    %   which is a standard forward IVP solved with ode45.
    %   A(t) and B(t) are obtained by central finite differences of f(x,u).
    %
    % ── Usage example ────────────────────────────────────────────────────────
    %   params.M = 1.0; params.m = 0.1; params.L = 0.5; params.g = 9.81;
    %   Q = diag([10, 1, 100, 1]);
    %   R = 0.01;
    %   tvlqr = computeTVLQR(t_nom, x_nom, u_nom, Q, R, params);
    %
    % ─────────────────────────────────────────────────────────────────────────
    
    %% ── 0. Input validation ──────────────────────────────────────────────────
    
    assert(isrow(t_nom),           't_nom must be a row vector [1 x N]');
    assert(size(x_nom,1) == 4,     'x_nom must be [4 x N]');
    assert(size(u_nom,1) == 1,     'u_nom must be [1 x N]');
    assert(size(x_nom,2) == length(t_nom), 'x_nom and t_nom must have same N');
    assert(size(u_nom,2) == length(t_nom), 'u_nom and t_nom must have same N');
    assert(issymmetric(Q) && all(eig(Q) > 0), 'Q must be symmetric positive definite');
    assert(isscalar(R)    && R > 0,            'R must be a positive scalar');
    
    nx = 4;   % state dimension
    nu = 1;   % input dimension
    N  = length(t_nom);
    t0 = t_nom(1);
    tf = t_nom(end);
    
    fprintf('[computeTVLQR] N=%d knots,  t in [%.3f, %.3f] s\n', N, t0, tf);
    
    %% ── 1. Build interpolants for x_nom(t) and u_nom(t) ─────────────────────
    % pchip: shape-preserving, no overshoot -- important for dynamics eval.
    
    x_interp = @(t) interp1(t_nom, x_nom', t, 'pchip', 'extrap')';  % [4 x 1]
    u_interp = @(t) interp1(t_nom, u_nom', t, 'pchip', 'extrap');    % scalar
    
    %% ── 2. Numerical linearisation helper ───────────────────────────────────
    % Central finite differences.  eps_fd = 1e-6 gives ~1e-12 truncation
    % error for smooth cart-pole dynamics.
    
    eps_fd = 1e-6;
    
        function [A, B] = linearise(t)
            xk = x_interp(t);
            uk = u_interp(t);
            A  = zeros(nx, nx);
            for i = 1:nx
                ei      = zeros(nx, 1); ei(i) = eps_fd;
                A(:,i)  = (cartPoleDynamics(xk+ei, uk, params) ...
                         - cartPoleDynamics(xk-ei, uk, params)) / (2*eps_fd);
            end
            B = (cartPoleDynamics(xk, uk+eps_fd, params) ...
               - cartPoleDynamics(xk, uk-eps_fd, params)) / (2*eps_fd);
        end
    
    %% ── 3. Solve the RDE backward via forward integration in tau = tf - t ───
    %
    %  The standard RDE in time t (from usual literature):
    %
    %  -dS/dt = A(t)'S + S A(t) - S B(t) R^{-1} B(t)' S + Q        ...(*)
    %
    %  ode45 only integrates forward, so we substitute tau = tf - t:
    %
    %    dS/dtau = dS/dt * dt/dtau = dS/dt * (-1) = -dS/dt
    %
    %  Negating (*) gives the forward-tau form -- LHS sign flips:
    %
    %   +dS/dtau = A(tf-tau)' S + S A(tf-tau) - S B R^{-1} B' S + Q  ...(**)
    %
    %  This is why there is NO leading minus sign in rde_rhs below.
    %  It is NOT a bug -- it is the result of the time-reversal substitution.
    %
    %  Initial condition at tau=0 (i.e. t=tf):  S(tf) = Q  (terminal cost).
    %  ode45 integrates (**) forward from tau=0 to tau=tf-t0.

    %
    %  State vector for ode45: s = vec(S), length 16.
    %  The RDE in forward-tau form:
    %    dS/dtau = A(tf-tau)' S + S A(tf-tau) - S B(tf-tau) R^{-1} B(tf-tau)' S + Q
    
    Rinv = 1 / R;
    
        function dSdtau = rde_rhs(tau, s_vec)
            t_curr    = tf - tau;                  % map back to physical time
            S_curr    = reshape(s_vec, nx, nx);    % 16-vec -> 4x4
            [A, B]    = linearise(t_curr);
            BRinvBt   = B * (Rinv * B');           % [4x4], outer product for nu=1
            dS        = A' * S_curr + S_curr * A - S_curr * BRinvBt * S_curr + Q;
            dSdtau    = dS(:);                     % 4x4 -> 16-vec
        end
    
    % Initial condition: S(tf) = Q  (terminal cost = Q, same as running cost)
    S0_vec = Q(:);
    
    % Integrate from tau=0 (t=tf) to tau=(tf-t0) (t=t0)
    tau_span = [0, tf - t0];
    
    ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, ...
                      'MaxStep', (tf-t0)/200);
    
    fprintf('[computeTVLQR] Integrating RDE backward in time ...\n');
    [tau_sol, S_sol] = ode45(@rde_rhs, tau_span, S0_vec, ode_opts);
    % tau_sol : [M x 1],  S_sol : [M x 16]
    
    %% ── 4. Map solution back to physical time and evaluate at t_nom ──────────
    %
    %  tau_sol is at ode45's adaptive steps, not at t_nom knots.
    %  We interpolate S(tau) -> S(t) by:
    %    tau = tf - t   =>   t = tf - tau
    %  tau_sol is increasing (0 -> tf-t0), so t_from_tau = tf - tau_sol
    %  is decreasing (tf -> t0). We flip both arrays so t increases.
    
    t_from_tau  = tf - tau_sol;          % [M x 1], decreasing
    t_ascending = flip(t_from_tau);      % [M x 1], increasing
    S_ascending = flip(S_sol, 1);        % [M x 16], matching flip
    
    % Interpolate onto original t_nom knots
    S_at_tnom = interp1(t_ascending, S_ascending, t_nom', 'pchip');
    % S_at_tnom : [N x 16]
    
    %% ── 5. Compute K(t) = R^{-1} B(t)' S(t) at each knot point ─────────────
    
    S_out = zeros(nx, nx, N);   % [4 x 4 x N]
    K_out = zeros(nu, nx, N);   % [1 x 4 x N]
    
    for k = 1:N
        Sk          = reshape(S_at_tnom(k,:), nx, nx);
        [~, Bk]     = linearise(t_nom(k));
        S_out(:,:,k) = Sk;
        K_out(:,:,k) = Rinv * (Bk' * Sk);    % [1x4]
    end
    
    %% ── 6. Sanity checks ─────────────────────────────────────────────────────
    
    % Check S(t) remains positive definite across all knots
    min_eigs = zeros(1, N);
    for k = 1:N
        min_eigs(k) = min(eig(S_out(:,:,k)));
    end
    
    if any(min_eigs <= 0)
        warning('[computeTVLQR] S(t) lost positive definiteness at %d knot(s). ', ...
                'Check Q, R, or trajectory validity. min eig = %.3e', ...
                sum(min_eigs <= 0), min(min_eigs));
    else
        fprintf('[computeTVLQR] S(t) is PD at all knots. min eig = %.4f\n', ...
                min(min_eigs));
    end
    
    % Check K is finite
    if any(~isfinite(K_out(:)))
        warning('[computeTVLQR] K(t) contains non-finite values.');
    end
    
    fprintf('[computeTVLQR] Done. Struct fields: t, S [4x4xN], K [1x4xN], Q, R\n');
    
    %% ── 7. Package output struct ─────────────────────────────────────────────
    
    tvlqr.t = t_nom;      % [1 x N]
    tvlqr.S = S_out;      % [4 x 4 x N]
    tvlqr.K = K_out;      % [1 x 4 x N]
    tvlqr.Q = Q;          % [4 x 4]
    tvlqr.R = R;          % scalar

end % ── end computeTVLQR ───────────────────────────────────────────────────

%% ═══════════════════════════════════════════════════════════════════════════
function f = cartPoleDynamics(x, u, p)
% CARTPOLEDYNAMICS  Evaluate f(x,u) for the cart-pole system.
%
%   x = [x_cart; v_x; theta; omega]
%   u = horizontal force (N), scalar
%   p = struct with fields: M, m, L, g

    v_x   = x(2);
    theta = x(3);
    omega = x(4);

    s = sin(theta);
    c = cos(theta);
    D = p.M + p.m * s^2;          % denominator, always > 0

    f = [
        v_x;
        (u  +  p.m*p.L*omega^2*s  +  p.m*p.g*s*c)            / D;
        omega;
        (-u*c  -  p.m*p.L*omega^2*s*c  -  (p.M+p.m)*p.g*s)   / (p.L * D)
    ];
end