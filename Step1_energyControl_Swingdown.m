clc; clearvars; close all;
params = init_params();

%% Compute nominal trajectory and associated nominal control
[t_nom, x_nom, u_nom] = generate_nominal_trajectory_and_input(params);

visualize_state_trajectory_and_input_history(t_nom, x_nom, u_nom, params);

t_nom = t_nom'; x_nom = x_nom'; u_nom = u_nom';
% control law function handle
control_law_fn_handle = @(t, x) energy_shaping_law(x, params);
save('./precomputedData/swing_down/nominal_trajectory_and_input.mat', ...
        "t_nom", "x_nom","u_nom","params","control_law_fn_handle");

keyboard

clc; clearvars; close all
analysis_mode = 'swing-down'; numSamples = 1000;
run('./utils/checkClosedLoop_MCRollouts.m');

% Add empirically-observed min-max bounds to the file
save('./precomputedData/swing_down/nominal_trajectory_and_input.mat', "minmaxBounds", "numSamples", '-append');

%% Function definitions

function params = init_params()
    params.M = 1.0; params.m = 0.1; params.L = 0.5; params.g = 9.81;
    params.F_max = 10;
    params.tspan = [0 15];
    params.num_knot_points = 100;

    % Initial and Final Centers
    params.x0 = [0; 0; pi; 0];
    params.xf = [0; 0; 0; 0];
    
    % LQR gains for attractor at top
    params.Q = diag([10, 1, 100, 1]);
    params.R = 5;

    %control law gains
    % params.gains.k_e = 2; params.gains.k_p = 1.0; params.gains.k_d = 0.5;
    params.gains.k_e = 1; params.gains.mu = 1;   
    params.gains.K = compute_attractor_gain(params);
    % params.initial_impulse = -0.01; %in N
    params.initial_perturbation = 1e-3; %to kickstart the swing-up controller
    params.controller_switch_threshold = 0.2*pi;
end

function x_dot = cartpole_dynamics(t, x, u, params)
    
    s = sin(x(3)); c = cos(x(3));
    v = x(2); omega = x(4);
    denom = params.M + params.m * s^2;
    
    f2 = (u + params.m*params.L*omega^2*s + params.m*params.g*s*c) / denom;
    f4 = (-u*c - params.m*params.L*omega^2*s*c - (params.M+params.m)*params.g*s) / (params.L * denom);
    x_dot = [v; f2; omega; f4];
end

function K = compute_attractor_gain(params)
    M = params.M; m= params.m; L = params.L; g = params.g;

    A = [0 1 0 0; 0 0 m*g/M 0; 0 0 0 1; 0 0 -(M+m)*g/(M*L) 0];
    B = [0; 1/M; 0; -1/(M*L)];

    Q = params.Q;   % penalise theta heavily
    R = params.R;

    K = lqr(A,B,Q,R);
end

function [t_nom, x_nom, u_nom] = generate_nominal_trajectory_and_input(params)

    tspan  = params.tspan;
    x_init = params.x0 + params.initial_perturbation * randn(size(params.x0));
    x_init(3) = wrapToPi(x_init(3));   % ensure initial theta is wrapped

    t_nom = []; x_nom = []; u_nom = [];
    t_start = tspan(1);
    t_end = tspan(2);
    x_curr   = x_init;

    opts = odeset('RelTol',1e-6,'AbsTol',1e-8, 'Events', @wrap_event);

    while t_start < t_end

        [t_seg, x_seg, te, xe, ~] = ode15s(@(t,x) cartpole_dynamics(t, x, ...
                       energy_shaping_law(x, params), params), ...
                            [t_start, t_end], x_curr, opts);

        % Reconstruct u for this segment
        u_seg = zeros(length(t_seg), 1);
        for i = 1:length(t_seg)
            u_seg(i) = energy_shaping_law(x_seg(i,:)', params);
        end

        % Append segment (skip duplicate point on restart)
        if isempty(t_nom)
            t_nom = t_seg;
            x_nom = x_seg;
            u_nom = u_seg;
        else
            t_nom = [t_nom; t_seg(2:end)];
            x_nom = [x_nom; x_seg(2:end,:)];
            u_nom = [u_nom; u_seg(2:end)];
        end

        % If no event fired, integration reached t_end -- done
        if isempty(te)
            break;
        end

        % Event fired: wrap theta and restart
        t_start = te(end);
        x_curr   = xe(end,:)';
        x_curr(3) = wrapToPi(x_curr(3));   % wrap theta to [-pi, pi]
    end

    %interpolate/extrapolate to match the number of knot points
    t_nom_temp = linspace(t_start,t_end,params.num_knot_points);
    x_nom = interp1(t_nom, x_nom, t_nom_temp, 'pchip', 'extrap');
    u_nom = interp1(t_nom, u_nom, t_nom_temp, 'pchip', 'extrap')';
    t_nom = t_nom_temp';
end

% ── Event: theta crosses +pi or -pi ──────────────────────────────
function [val, isterminal, direction] = wrap_event(t, x)
    % Fires when theta - pi = 0 (crossing +pi)
    % or     when theta + pi = 0 (crossing -pi)
    val        = [x(3) - pi;    % crossing +pi
                  x(3) + pi];   % crossing -pi
    isterminal = [1; 1];        % stop integration at either crossing
    direction  = [0; 0];        % fire on any crossing direction
end


function u = energy_shaping_law(x, params)
    
    M = params.M; m= params.m; L = params.L; g = params.g;
    theta = wrapToPi(x(3)); omega = x(4);

    % % A constant initial impulse to kickstart the energy pumping-based swing-up strategy
    % if norm(x - params.x0) < 1e-3
    %     u = params.initial_impulse;
    %     return
    % end

    % Energy of pendulum: E = 0.5*m*L^2*w^2 - m*g*L*cos(theta)
    E = 0.5*m*L^2*omega^2 - m*g*L*cos(theta);
    E_down = -m*g*L;
    
    k_e = params.gains.k_e; mu = params.gains.mu;
    K = params.gains.K;

    % Pump energy based on states
    if abs(theta - params.xf(3)) > params.controller_switch_threshold
        % u = k_e * (E - E_down) * omega * cos(theta);
        % u = k_e * (M + m*sin(theta)^2) * (E - E_down)*omega*cos(theta);
        % u = k_e * (M + m*sin(theta)^2) * ((E - E_down)*omega*cos(theta) - mu*x(2)) - ...
        %     m*sin(theta) * (L*omega^2 + g*cos(theta)); %exact control law considering feedback terms as well
        u = k_e * ((E - E_down)*omega*cos(theta) - mu*x(2)) - ...
            m*sin(theta) * (L*omega^2 + g*cos(theta)); %exact control law considering feedback terms as well
    else
        u = K*(params.xf - x);
    end

    %Apply input saturations
    u = max(min(u, params.F_max), -params.F_max);
end


function visualize_state_trajectory_and_input_history(t_nom, x_nom, u_nom, params)

    %state trajectory history
    figure;
    subplot(2,2,1); hold on; grid on;
    plot(t_nom, x_nom(:,1), 'b', 'LineWidth', 1.5);
    xlabel('t (s)'); ylabel('Cart Position (m)');
    
    subplot(2,2,2); hold on; grid on;
    plot(t_nom, x_nom(:,2), 'b', 'LineWidth', 1.5);
    xlabel('t (s)'); ylabel('Cart Velocity (m/s)');
    
    subplot(2,2,3); hold on; grid on;
    %x_nom(:,3) = wrapToPi(x_nom(:,3));
    plot(t_nom, rad2deg(x_nom(:,3)), 'b', 'LineWidth', 1.5);
    yline(0, 'k--', 'Hanging-down');
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

    % % Lyapunov function and Energy over time
    % % Energy of pendulum: E = 0.5*m*L^2*w^2 - m*g*L*cos(theta)
    % M = params.M; m= params.m; L = params.L; g = params.g;
    % E = 0.5*m*L^2*x_nom(:,4).^2 - m*g*L*cos(x_nom(:,3));
    % E_down = -m*g*L;
    % 
    % V = 0.5*(E-E_down).^2 + 0.5*params.gains.mu*x_nom(:,2).^2;
    % 
    % %Lyapunov function and energy error
    % figure;
    % subplot(2,1,1); hold on; grid on;
    % plot(t_nom, (E-E_down), 'k', 'LineWidth', 1.5);
    % xlabel('t (s)'); ylabel('Energy difference (J)');
    % 
    % subplot(2,1,2); hold on; grid on;
    % plot(t_nom, V, 'k', 'LineWidth', 1.5);
    % xlabel('t (s)'); ylabel('Lyapunov function value, V');

    % input profile
    figure; grid on; hold on
    plot(t_nom, u_nom, 'k-.', 'LineWidth', 1.75);
    ylim([-8, 8]);
    xlabel('t (s)'); ylabel('F (N)');
    title('Input profile');

end