clc; clearvars; close all;

params = init_params();

%% Compute nominal trajectory and associated nominal control
[t_nom, x_nom, u_nom] = generate_nominal_trajectory_and_input(params);

visualize_state_trajectory_and_input_history(t_nom, x_nom, u_nom);

t_nom = t_nom'; x_nom = x_nom'; u_nom = u_nom';
save('./precomputedData/swing_up/nominal_trajectory_and_input.mat',"t_nom", "x_nom","u_nom","params");

%% Compute a stabilizing TVLQR feedback controller
% Load the nominal trajectory and feedforward control, and 
load('./precomputedData/swing_up/nominal_trajectory_and_input.mat');

% Run and save only once -- after that only use stored values 
% Compute the TVLQR gains
Q = diag([10, 1, 100, 1]);   % penalise theta heavily
R = 0.01;

tvlqr = computeTVLQR(t_nom, x_nom, u_nom, Q, R, params);
save('./precomputedData/swing_up/TVLQR_gains_and_cost_matrices.mat', "tvlqr")

%% Function definitions

function params = init_params()
    params.M = 1.0; params.m = 0.1; params.L = 0.5; params.g = 9.81;
    params.F_max = 10;
    params.tspan = [0 15]; % 6 seconds for a "gentle" swing up
    % Initial and Final Centers
    params.x0 = [0; 0; 0; 0];
    params.xf = [0; 0; pi; 0];
    % Set definitions (Ellipsoid S-matrices)
    params.P_0 = 100*diag([10, 0.5, 10, 0.25]); % Example X0 size
    params.P_f = diag([1, 1, 0.5, 0.5]); % Requirement for Xf

    %control law gains
    % params.gains.k_e = 2; params.gains.k_p = 1.0; params.gains.k_d = 0.5;
    params.gains.k_e = 2; params.gains.k_p = 0; params.gains.k_d = 0;
    params.gains.K = 0.5*[-10.0000  -16.2819   91.7720   22.6933];
    params.initial_impulse = 3; %in N
    params.controller_switch_angle = 0.8*pi;
end

function dx = cartpole_dynamics(t, x, u, params)
    s = sin(x(3)); c = cos(x(3));
    v = x(2); omega = x(4);
    denom = params.M + params.m * s^2;
    
    f2 = (u + params.m*params.L*omega^2*s + params.m*params.g*s*c) / denom;
    f4 = (-u*c - params.m*params.L*omega^2*s*c - (params.M+params.m)*params.g*s) / (params.L * denom);
    dx = [v; f2; omega; f4];
end

function [t_nom, x_nom, u_nom] = generate_nominal_trajectory_and_input(params)
    tspan = params.tspan;
    x_init = params.x0;

    % Energy Shaping Controller
    ctrl = @(t, x) energy_shaping_law(x, params);
    
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t_nom, x_nom] = ode45(@(t, x) cartpole_dynamics(t, x, ctrl(t, x), params), tspan, x_init, options);
    
    % Reconstruct u
    u_nom = zeros(length(t_nom), 1);
    for i = 1:length(t_nom)
        u_nom(i) = ctrl(t_nom(i), x_nom(i,:)');
    end
end

function u = energy_shaping_law(x, params)
    
    % A constant initial impulse to kickstart the energy pumping-based swing-up strategy
    if norm(x - params.x0) < 1e-3
        u = params.initial_impulse;
        return
    end

    % Energy of pendulum: E = 0.5*m*L^2*w^2 - m*g*L*cos(theta)
    E = 0.5*params.m*params.L^2*x(4)^2 - params.m*params.g*params.L*cos(x(3));
    E_up = params.m*params.g*params.L;
    
    k_e = params.gains.k_e; k_p = params.gains.k_p; k_d = params.gains.k_d;
    K = params.gains.K;

    % Pump energy based on velocity and position
    if abs(x(3)) < params.controller_switch_angle
        u = k_e * (E - E_up) * x(4) * cos(x(3)); 
    else
        u = K*(params.xf - x);
    end

    % PD to keep cart near origin: not activated currently! (zero-gains)
    u = u - k_p*x(1) - k_d*x(2);
    u = max(min(u, params.F_max), -params.F_max);
end

function [A_func, B_func] = get_symbolic_jacobians()
    % Define symbolic variables
    syms x v theta omega F M m L g real
    
    % Dynamics (Your provided equations)
    s = sin(theta); c = cos(theta);
    denom = M + m*s^2;
    
    f1 = v;
    f2 = (F + m*L*omega^2*s + m*g*s*c) / denom;
    f3 = omega;
    f4 = (-F*c - m*L*omega^2*s*c - (M + m)*g*s) / (L * denom);
    
    f_sym = [f1; f2; f3; f4];
    state = [x; v; theta; omega];
    
    % Compute Jacobians
    A_sym = jacobian(f_sym, state);
    B_sym = jacobian(f_sym, F);
    
    % Convert to fast MATLAB functions
    % We pass params as a structure to the function
    A_func = matlabFunction(A_sym, 'Vars', {x, v, theta, omega, F, M, m, L, g});
    B_func = matlabFunction(B_sym, 'Vars', {x, v, theta, omega, F, M, m, L, g});
end

% %Changed this into symbolic for numeric stability -- so not required
% function [A, B] = get_jacobians(x, u, params)
%     % Finite difference or symbolic derivatives of cartpole_dynamics
%     eps = 1e-6; nx = 4; nu = 1;
%     A = zeros(nx, nx); B = zeros(nx, nu);
%     for i = 1:nx
%         x_plus = x; x_plus(i) = x_plus(i) + eps;
%         A(:,i) = (cartpole_dynamics(0, x_plus, u, params) - cartpole_dynamics(0, x, u, params))/eps;
%     end
%     u_plus = u + eps;
%     B = (cartpole_dynamics(0, x, u_plus, params) - cartpole_dynamics(0, x, u, params))/eps;
% end

function K_list = compute_tvlqr_gains(t_nom, x_nom, u_nom, params)
    Q = diag([10, 1, 50, 1]); R = 10; Pf = params.P_f;
    N = length(t_nom);
    
    % Get symbolic jacobian functions
    [A_func, B_func] = get_symbolic_jacobians();

    % For simplicity in this pass, we solve a sequence of discrete LQR
    % In a stiff system, you'd solve the Differential Riccati Equation.
    K_list = zeros(N, 4);
    for i = N:-1:1
        % [A, B] = get_jacobians(x_nom(i,:)', u_nom(i), params);
        A = A_func(x_nom(i,1), x_nom(i,2), x_nom(i,3), x_nom(i,4), u_nom(i),...
                        params.M, params.m, params.L, params.g);
        B = B_func(x_nom(i,1), x_nom(i,2), x_nom(i,3), x_nom(i,4), u_nom(i),...
                        params.M, params.m, params.L, params.g);

        % This is a quasi-static approximation for the gains
        if i == N %use P_f for terminal time only
            [K_list(i,:), ~, ~] = lqr(A, B, Pf, R); 
        else %otherwise use Q stage cost matrix
            [K_list(i,:), ~, ~] = lqr(A, B, Q, R); 
        end
    end
end

%To Do: enforce this only at level-set boundary: X^T S X = 1
function [t_S, S_history, S_vec_hist] = propagate_reachability(t_nom, x_nom, u_nom, K_feedback, params)
    
    n_x = size(x_nom,2);
    S0_vec = params.P_0(:);
    
    % Get symbolic jacobians
    [A_func, B_func] = get_symbolic_jacobians();

    % Use ode15s for stiffness handling
    options = odeset('RelTol', 1e-6);
    [t_S, S_vec_hist] = ode15s(@(t, s) lyapunov_rhs(t, s, t_nom, x_nom, u_nom, A_func, B_func, K_feedback, params), ...
                                [t_nom(1) t_nom(end)], S0_vec, options);
    
    N = size(S_vec_hist,1);
    S_history = nan(n_x,n_x,N);
    for k = 1:N
        S_history(:,:,k) = reshape(S_vec_hist(k, :), [n_x, n_x]);
    end

end

function ds_vec = lyapunov_rhs(t, s_vec, t_nom, x_nom, u_nom, A_func, B_func, K_feedback, params)
    % Interpolate nominal data
    xn = interp1(t_nom, x_nom, t, "pchip")';
    un = interp1(t_nom, u_nom, t, "pchip");
    K = interp1(t_nom, K_feedback, t, "pchip");
    
    S = reshape(s_vec, [4, 4]);
    % [A, B] = get_jacobians(xn, un, params);
    A = A_func(xn(1), xn(2), xn(3), xn(4), un, params.M, params.m, params.L, params.g);
    B = B_func(xn(1), xn(2), xn(3), xn(4), un, params.M, params.m, params.L, params.g);
    
    A_cl = A - B * K;
    
    % S-propagation: dS = -(A_cl'*S + S*A_cl)
    dS = -(A_cl' * S + S * A_cl);

    ds_vec = dS(:);
end

function visualize_state_trajectory_and_input_history(t_nom, x_nom, u_nom)
    
    %state trajectory history
    figure;
    subplot(2,2,1); hold on; grid on;
    plot(t_nom, x_nom(:,1), 'b', 'LineWidth', 1.5);
    xlabel('t (s)'); ylabel('Cart Position (m)');
    
    subplot(2,2,2); hold on; grid on;
    plot(t_nom, x_nom(:,2), 'b', 'LineWidth', 1.5);
    xlabel('t (s)'); ylabel('Cart Velocity (m/s)');
    
    subplot(2,2,3); hold on; grid on;
    plot(t_nom, rad2deg(x_nom(:,3)), 'b', 'LineWidth', 1.5);
    yline(180, 'k--', 'Upright');
    xlabel('t (s)'); ylabel('\theta (rad)');
    
    subplot(2,2,4); hold on; grid on;
    plot(t_nom, x_nom(:,4), 'b', 'LineWidth', 1.5);
    xlabel('t (s)'); ylabel('\theta dot (rad/s)');
    
    sgtitle('Cart-Pole State Trajectories');
    
    % x-theta trajectory
    figure; grid on; hold on; axis equal;
    plot(x_nom(:,1), x_nom(:,3), 'b--', 'LineWidth', 1.5);
    xlabel('Cart Position (m)'); ylabel('\theta (rad)');
    plot(x_nom(1,1), x_nom(1,3), 'sg', 'MarkerSize', 7, 'LineWidth', 1.5);
    plot(x_nom(end,1), x_nom(end,3), 'xr', 'MarkerSize', 7, 'LineWidth', 1.5);
    title('x-\theta');

    % input profile
    figure; grid on; hold on
    plot(t_nom, u_nom, 'k-.', 'LineWidth', 1.75);
    xlabel('t (s)'); ylabel('F (N)');
    title('Input profile');

end

function visualize_reachability(t_nom, x_nom, t_S, S_history)
    
    N   = size(x_nom,1); %number of knot points (discretization)
    n_x = size(x_nom,2); %state-space dimensionality

    % C selects state 1 (x) and state 3 (theta)
    C = [1, 0, 0, 0; 
         0, 0, 1, 0]; 
    
    figure('Color', 'w', 'Name', 'Reachability Analysis: x-theta Projection');
    hold on; grid on;
    
    % 1. Plot the Nominal Trajectory in x-theta space
    plot(x_nom(:,1), x_nom(:,3), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Nominal Path');
    
    % 2. Plot Ellipsoids at discrete intervals
    num_ellipsoids = round(0.1*N);
    indices = round(linspace(1, size(S_history, 1), num_ellipsoids));
    
    colors = winter(num_ellipsoids);
    
    for i = 1:num_ellipsoids
        idx = indices(i);
        t_curr = t_S(idx);
        
        % Reshape and invert to get the "Dual" (Covariance) matrix
        S_4d = reshape(S_history(idx, :), [n_x, n_x]);
        
        % Numerical check: ensure S is positive definite before inverting
        [V, D] = eig(S_4d);
        D_clean = max(D, 1e-6); % Regularize for inversion
        S_clean = V * D_clean * V';
        P_4d = inv(S_clean);
        
        % Project onto 2D
        P_2d = C * P_4d * C';
        S_2d = inv(P_2d);
        
        % Get the center at this time (interpolate from x_nom)
        center = interp1(t_nom, x_nom, t_curr, "pchip")';
        c_2d = C * center;
        
        % Draw the 2D Ellipse
        draw_ellipse(c_2d, S_2d, colors(i, :), 0.3);
    end
    
    % Formatting
    xlabel('Cart Position x (m)');
    ylabel('Pendulum Angle \theta (rad)');
    title('Projected Reachable Funnel (x-\theta plane)');
    % legend('Location', 'best');
    axis equal;
end

function draw_ellipse(center, S, color, alpha)
    % Generate points for a unit circle
    theta = linspace(0, 2*pi, 100);
    unit_circle = [cos(theta); sin(theta)];
    
    % Transform circle to ellipsoid: x = center + L*unit_circle
    % where L is the Cholesky factor of the inverse shape matrix (P)
    P = inv(S);
    try
        L = chol(P, 'lower');
        points = center + L * unit_circle;
        fill(points(1,:), points(2,:), color, 'FaceAlpha', alpha, ...
             'EdgeColor', color, 'HandleVisibility', 'off');
    catch
        % If not positive definite, plot a marker instead
        plot(center(1), center(2), 'rx');
    end
end