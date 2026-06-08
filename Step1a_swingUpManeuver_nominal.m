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

run('./utils/checkClosedLoop_MCRollouts.m');

%% Function definitions

function params = init_params()
    params.M = 1.0; params.m = 0.1; params.L = 0.5; params.g = 9.81;
    params.F_max = 10;
    params.tspan = [0 15]; % 6 seconds for a "gentle" swing up
    % Initial and Final Centers
    params.x0 = [0; 0; 0; 0];
    params.xf = [0; 0; pi; 0];
    
    % LQR gains
    params.Q = diag([10, 1, 100, 1]);
    params.R = 0.5;

    %control law gains
    % params.gains.k_e = 2; params.gains.k_p = 1.0; params.gains.k_d = 0.5;
    params.gains.k_e = 2; params.gains.k_p = 0; params.gains.k_d = 0;
    % params.gains.K = 0.5*[-10.0000  -16.2819   91.7720   22.6933];
    params.gains.K = compute_attractor_gain(params);
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

function K = compute_attractor_gain(params)
    M = params.M; m= params.m; L = params.L; g = params.g;

    A = [0 1 0 0; 0 0 m*g/M 0; 0 0 0 1; 0 0 (M+m)*g/(M*L) 0];
    B = [0; 1/M; 0; 1/(M*L)];

    Q = params.Q;   % penalise theta heavily
    R = params.R;

    K = lqr(A,B,Q,R);
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