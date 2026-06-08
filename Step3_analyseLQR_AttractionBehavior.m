clc; clearvars; close all;
addpath('./lib/');

%% Analyze LQR attraction behavior
params = init_params('bottom'); %OPTIONS: 'top' or 'bottom'
x_init = sample_points_from_ellipsoid(params.initial_set_matrix, params.xf, 1, 'boundary'); %OPTIONS: 'interior' or 'boundary'
[t_nom, x_nom, u_nom] = analyze_trajectory_and_input(params, x_init);

visualize_state_trajectory_and_input_history(t_nom, x_nom, u_nom, params);

function [t_nom, x_nom, u_nom] = analyze_trajectory_and_input(params, x_init)
    tspan = params.tspan;

    % Energy Shaping Controller
    ctrl = @(t, x) LQR_attractor_at_top(x, params);
    
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t_nom, x_nom] = ode45(@(t, x) cartpole_dynamics(t, x, ctrl(t, x), params), tspan, x_init, options);
    
    % Reconstruct u
    u_nom = zeros(length(t_nom), 1);
    for i = 1:length(t_nom)
        u_nom(i) = ctrl(t_nom(i), x_nom(i,:)');
    end
end

function u = LQR_attractor_at_top(x, params)

    K = params.gains.K;
    u = K*(params.xf - x);
    u = max(min(u, params.F_max), -params.F_max);
end

function params = init_params(mode)
    params.M = 1.0; params.m = 0.1; params.L = 0.5; params.g = 9.81;
    params.F_max = 10;
    params.tspan = [0 15]; % 6 seconds for a "gentle" swing up
    
    if strcmpi(mode, 'top')
        % At the top (theta = pi)
        params.xf = [0; 0; pi; 0]; %final state
    elseif strcmpi(mode, 'bottom')
        % At the bottom (theta = 0)
        params.xf = [0; 0; 0; 0]; %final state
    else
        error('Analysis mode must be either top or bottom');
    end
    params.mode = mode;    
    
    % LQR gains for attractor at top
    params.Q = diag([10, 1, 100, 1]); %penalize theta heavily
    params.R = 0.5;

    %control law gains
    % params.gains.K = 0.5*[-10.0000  -16.2819   91.7720   22.6933];
    [params.gains.K, params.S] = compute_attractor_gain(params);

    params.initial_set_matrix = 0.1*params.S; %some scaling of Riccati solution

    % params.initial_impulse = 3; %in N
    % params.controller_switch_angle = 0.8*pi;
end

function dx = cartpole_dynamics(t, x, u, params)
    s = sin(x(3)); c = cos(x(3));
    v = x(2); omega = x(4);
    denom = params.M + params.m * s^2;
    
    f2 = (u + params.m*params.L*omega^2*s + params.m*params.g*s*c) / denom;
    f4 = (-u*c - params.m*params.L*omega^2*s*c - (params.M+params.m)*params.g*s) / (params.L * denom);
    dx = [v; f2; omega; f4];
end

function [K, S] = compute_attractor_gain(params)
    M = params.M; m= params.m; L = params.L; g = params.g;

    if strcmpi(params.mode, 'top')
        % At the top (theta = pi)
        A = [0 1 0 0; 0 0 m*g/M 0; 0 0 0 1; 0 0 (M+m)*g/(M*L) 0];
        B = [0; 1/M; 0; 1/(M*L)];
    elseif strcmpi(params.mode, 'bottom')
        % At the bottom (theta = 0)
        A = [0 1 0 0; 0 0 m*g/M 0; 0 0 0 1; 0 0 -(M+m)*g/(M*L) 0];
        B = [0; 1/M; 0; -1/(M*L)];
    else
        error('Analysis mode must be either top or bottom');
    end

    Q = params.Q;
    R = params.R;

    [K, S, ~] = lqr(A,B,Q,R);
end

function points = sample_points_from_ellipsoid(M, xc, N, type)
% SAMPLE_ELLIPSE_GENERAL Samples N points from an n-dimensional ellipse
% Define by: (x - xc)' * M * (x - xc) <= 1
%
% Inputs:
%   M    - n x n symmetric positive-definite matrix
%   xc   - n x 1 column vector representing the center of the ellipse
%   N    - Number of points to sample (scalar integer)
%   type - String, either 'interior' or 'boundary'
%
% Output:
%   points - N x n matrix where each row is an n-dimensional sampled point

    % 1. Validate inputs and dimensions
    [n, m] = size(M);
    if n ~= m || any(eig(M) <= 0)
        error('M must be a square, symmetric positive-definite matrix.');
    end
    
    xc = xc(:); % Ensure xc is a column vector
    if length(xc) ~= n
        error('Dimensions of M and xc must match.');
    end
    
    type = lower(type);
    if ~strcmp(type, 'interior') && ~strcmp(type, 'boundary')
        error('Type must be either ''interior'' or ''boundary''.');
    end

    % 2. Compute Cholesky decomposition of the inverse matrix
    % M^-1 = L * L' -> L maps a unit hypersphere to the target hyper-ellipse
    Minv = inv(M);
    L = chol(Minv, 'lower');

    % 3. Generate random points on an n-dimensional unit hypersphere surface
    % Standard normal distributions yield uniformly distributed directions
    z = randn(N, n); 
    norms = sqrt(sum(z.^2, 2));
    u_surface = z ./ norms; % Project points onto the exact surface (norm = 1)

    % 4. Apply radial scaling based on selection type
    if strcmp(type, 'interior')
        % In n-dimensions, volume scales with r^n.
        % To keep density uniform, we take the n-th root of a uniform variable.
        r = rand(N, 1).^(1 / n);
        u = u_surface .* r; % Scale points into the interior ball
    else
        u = u_surface; % Keep points on the exact boundary sphere
    end

    % 5. Transform unit ball/sphere points to the final hyper-ellipse
    % Transposed math: points = (L * u')' + xc' -> points = u * L' + xc'
    points = u * L' + xc';
end

function visualize_state_trajectory_and_input_history(t_nom, x_nom, u_nom, params)
    
    P = plottingFnsClass();
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
    % yline(180, 'k--', 'Upright');
    xlabel('t (s)'); ylabel('\theta (deg)');
    
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

    %plot initial sampling set
    projectionDims = [1 3];
    projected_ellipse_center = [params.xf(projectionDims(1)), params.xf(projectionDims(2))]';
    projected_inlet = P.project_ellipsoid_matrix_2D(params.initial_set_matrix, projectionDims); % Extract 2D covariance
    [eig_vec, eig_val] = eig(projected_inlet);
    
    theta = linspace(0, 2*pi, 100);
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(projected_ellipse_center(1) + rotated_ellipse(1, :), ...
         projected_ellipse_center(2) + rotated_ellipse(2, :), ...
         'k-.', 'LineWidth', 1.5);

    % input profile
    figure; grid on; hold on
    plot(t_nom, u_nom, 'k-.', 'LineWidth', 1.75);
    xlabel('t (s)'); ylabel('F (N)');
    title('Input profile');
end