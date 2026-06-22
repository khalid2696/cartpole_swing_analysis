% clc; clearvars; close all

%% Add directories
addpath('../lib/');

if ~exist('analysis_mode','var')
    analysis_mode = 'swing-down';
end

%% Load the nominal trajectory and control law
if strcmpi(analysis_mode, 'swing-up')
    load('../precomputedData/swing_up/nominal_trajectory_and_input.mat');
elseif strcmpi(analysis_mode, 'swing-down')
    load('../precomputedData/swing_down/nominal_trajectory_and_input.mat');
else
    error('Unsupported analysis mode')
end

fprintf('Running Monte Carlo rollouts for empirical analysis of the closed loop system.\nHang on..\n\n');

N = length(t_nom);

% -- Status message: Quantities at our disposal now -- %

% t_nom                 : time horizon (sampled)        : 1 x N
% x_nom                 : nominal state trajectory      : n_x x N
% control_law_fn_handle : parametrized control law dependent on state (x) only
% dynamicsFnHandle      : the function handle of the system dynamics

%% Specify parameters or Inherit them if they exist in the wrapper file

if ~exist('numSamples','var')
    numSamples = 1000; %default number of rollouts
end

%Sampling initial states from an initial ellipsoidal set
if strcmpi(analysis_mode, 'swing-up')
    initialStateSetMatrix = 121.4437 * eye(4); %based on funnel-outlets (at bottom) in the library
    finalStateSetCenter = [0 0 pi 0]';
    finalStateSetMatrix = 192.2777 * eye(4); %based on funnel-inlets (at top) in the library
elseif strcmpi(analysis_mode, 'swing-down')
    initialStateSetMatrix = 553.8465 * eye(4); %based on funnel-outlets (at top) in the library
    finalStateSetCenter = [0 0 0 0]';
    finalStateSetMatrix = 33.1433 * eye(4); %based on funnel-inlets (at bottom) in the library
else
    error('Unsupported analysis mode')
end
initialStateSetCenter = params.x0;

%% Monte Carlo rollouts
trajectories = cell(numSamples, 1);
inputProfiles = cell(numSamples, 1);

x_min_max = NaN(numSamples,2); %min-max bounds of x (to be used in motion planning code)
errors = zeros(numSamples, length(t_nom));

success = NaN(numSamples,1);

rollout_count = 1;
while rollout_count <= numSamples
    rollout_count

    %sample initial states at random
    x0 = sample_points_from_ellipsoid(initialStateSetMatrix, initialStateSetCenter, 1, 'boundary'); %Two options: 'interior' and 'boundary'
    
    [x_traj, control_input, x_bounds, error_norm, errorFlag] = forward_propagate(params, x0, ...
                                                    control_law_fn_handle, t_nom, x_nom, analysis_mode);
    if errorFlag 
        continue %repeat the MC rollout -- observed only in the swing-down maneuver (because of numerical ode error and integration error)
    end
    trajectories{rollout_count} = x_traj;
    inputProfiles{rollout_count} = control_input;
    x_min_max(rollout_count, :) = x_bounds;
    errors(rollout_count, :) = error_norm;

    % Compute the success of being within the user-specified terminal set
    success(rollout_count) = isContained(x_traj, finalStateSetCenter, finalStateSetMatrix);
    
    rollout_count = rollout_count+1;
end

successRate = mean(success);
disp('-- End of Monte Carlo forward rollouts --');
disp(' ');

%% Visualization 
disp('Plotting trajectories, input profiles, and metrics from MC rollouts..');
disp(' ');

%plot state trajectories
projectionDims = [1 3]; 
plot_xy_trajectories(trajectories, x_nom, initialStateSetMatrix, projectionDims);
plot_terminal_set(finalStateSetCenter, finalStateSetMatrix, projectionDims);
title(compose('Monte Carlo Rollouts, Success rate = %.2f from %d samples', successRate, numSamples));

%Plot input profile
plot_input_profiles(inputProfiles, u_nom, t_nom);
plot_error_metrics(errors, t_nom);

fprintf("Success rate of final state being within the specified terminal set is %.2f (from %d MC rollouts)\n",successRate,numSamples);

minmaxBounds = [min(x_min_max(:,1)), max(x_min_max(:,2))];
fprintf("\nMinimum x from %d MC rollouts: %0.3f", numSamples, minmaxBounds(1));
fprintf("\nMaximum x from %d MC rollouts: %0.3f", numSamples, minmaxBounds(2));
disp(' ');

%% Function defintions

function x_dot = cartpole_dynamics(t, x, u, params)
    s = sin(x(3)); c = cos(x(3));
    v = x(2); omega = x(4);
    denom = params.M + params.m * s^2;
    
    f2 = (u + params.m*params.L*omega^2*s + params.m*params.g*s*c) / denom;
    f4 = (-u*c - params.m*params.L*omega^2*s*c - (params.M+params.m)*params.g*s) / (params.L * denom);
    x_dot = [v; f2; omega; f4];
end

function [x_traj, control_input, x_bounds, errorNorms, errorFlag] = forward_propagate(params, x_init, ctrl, t_nom, x_nom, analysis_mode)
    
    errorFlag = 0;
    % swing-up and swing-down have some subtle ode implementation differences 
    % because of theta wrapping (0 to 2pi vs -pi to pi)
    switch analysis_mode
        case 'swing-up'
            options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); %'MaxStep',0.001
            [~, x_traj] = ode45(@(t, x) cartpole_dynamics(t, x, ctrl(t, x), params), t_nom, x_init, options);
            x_traj = x_traj';
        
            % Reconstruct u
            control_input = zeros(1, length(t_nom));
            errorNorms = zeros(length(t_nom),1);
            for k = 1:length(t_nom)
                control_input(k) = ctrl(t_nom(k), x_traj(:,k));
                errorNorms(k) = norm(x_traj(:, k) - x_nom(:, k));
            end
        case 'swing-down'
            x_init(3) = wrapToPi(x_init(3));   % ensure initial theta is wrapped
        
            time_instances = []; x_traj = []; control_input = [];
            t_start = t_nom(1);
            t_end = t_nom(end);
            x_curr   = x_init;
        
            opts = odeset('RelTol',1e-4,'AbsTol',1e-6, 'Events', @wrap_event);
            iter = 0;
            MAX_ITER = 100; % more than enough wrap events for any reasonable trajectory

            while t_start < t_end
                iter = iter + 1;
                if iter > MAX_ITER
                    warning('MC rollout %d: max wrap iterations reached at t=%.3f, breaking.', i, t_start);
                    errorFlag = 1;
                    %assign dummy values to the other outputs (vacuous)
                    x_traj = NaN; control_input = NaN;
                    x_bounds = [NaN, NaN]; errorNorms = NaN;
                    return
                end
                [t_seg, x_seg, te, xe, ~] = ode45(@(t,x) cartpole_dynamics(t, x, ...
                               ctrl(t, x), params), ...
                                    [t_start, t_end], x_curr, opts);
        
                % Reconstruct u for this segment
                u_seg = zeros(length(t_seg), 1);
                for i = 1:length(t_seg)
                    u_seg(i) = ctrl(t_seg(i), x_seg(i,:)');
                end
        
                % Append segment (skip duplicate point on restart)
                if isempty(time_instances)
                    time_instances = t_seg;
                    x_traj = x_seg;
                    control_input = u_seg;
                else
                    time_instances = [time_instances; t_seg(2:end)];
                    x_traj = [x_traj; x_seg(2:end,:)];
                    control_input = [control_input; u_seg(2:end)];
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
        
            % interpolate/extrapolate to match the number of knot points
            x_traj = interp1(time_instances, x_traj, t_nom, 'pchip', 'extrap');
            control_input = interp1(time_instances, control_input, t_nom, 'pchip', 'extrap')';

            % applying transpose to match convention
            x_traj = x_traj'; control_input = control_input'; 
            % Compute error norms
            errorNorms = sqrt(sum((x_traj-x_nom).^2, 1))';
            % errorNorms = zeros(length(t_nom),1);
            % for k = 1:length(t_nom)
            %     errorNorms(k) = norm(x_traj(:, k) - x_nom(:, k));
            % end
    end

    x_bounds = [min(x_traj(1,:)), max(x_traj(1,:))];
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

function points = sample_points_from_ellipsoid(M, xc, N, type)
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

% Outputs 1 if final state of rollout trajectory is within the user-specified terminal set
function success = isContained(x_traj, finalStateSetCenter, finalStateSetMatrix)
    if (x_traj(:,end) - finalStateSetCenter)' * finalStateSetMatrix * (x_traj(:,end) - finalStateSetCenter) < 1
        success = 1;
    else
        success = 0;
    end

end

function plot_xy_trajectories(trajectories, x_nom, initialSet, projectionDims)
    % Plot all trajectories and nominal trajectory
    figure; hold on; grid on; axis equal;
    
    P = plottingFnsClass();

    if nargin < 5
        projectionDims = [1 3];
    end
 
    plot(x_nom(projectionDims(1), :), x_nom(projectionDims(2), :), 'k--', 'LineWidth', 2);

    for i = 1:length(trajectories)
        plot(trajectories{i}(projectionDims(1), :), trajectories{i}(projectionDims(2), :), 'b-', 'LineWidth', 0.5);
        plot(trajectories{i}(projectionDims(1), 1), trajectories{i}(projectionDims(2), 1), 'sg');
    end
    
    %plot initial sampling set
    projected_ellipse_center = [x_nom(projectionDims(1),1), x_nom(projectionDims(2),1)]';
    projected_inlet = P.project_ellipsoid_matrix_2D(initialSet, projectionDims); % Extract 2D covariance
    [eig_vec, eig_val] = eig(projected_inlet);
    
    theta = linspace(0, 2*pi, 100);
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(projected_ellipse_center(1) + rotated_ellipse(1, :), ...
         projected_ellipse_center(2) + rotated_ellipse(2, :), ...
         'm-.', 'LineWidth', 1.5);

    xlabel(['x_{', num2str(projectionDims(1)), '}'])
    ylabel(['x_{', num2str(projectionDims(2)), '}'])

    legend('Nominal Trajectory','Rollout Trajectories','Sampled Initial states','Location','best');
    %hold off;
end

function plot_terminal_set(center, ellipsoidMatrix, projectionDims)
 
    P = plottingFnsClass();

    projected_ellipse_center = [center(projectionDims(1)), center(projectionDims(2))]';
    projected_matrix = P.project_ellipsoid_matrix_2D(ellipsoidMatrix, projectionDims); % Extract 2D covariance
    [eig_vec, eig_val] = eig(projected_matrix);
    
    theta = linspace(0, 2*pi, 100);
    ellipse_boundary = eig_val^(-1/2) * [cos(theta); sin(theta)];
    rotated_ellipse = eig_vec * ellipse_boundary;
    
    plot(projected_ellipse_center(1) + rotated_ellipse(1, :), ...
         projected_ellipse_center(2) + rotated_ellipse(2, :), ...
         'r-', 'LineWidth', 1.5, 'DisplayName', 'Terminal Set');
    plot(projected_ellipse_center(1), projected_ellipse_center(2), 'xk','MarkerSize',5, 'DisplayName', 'Terminal State');
end

function plot_state_trajectories(trajectories, t_nom, x_nom, stateDims)
    
    if nargin < 4
        stateDims = [1 2 3];
    end

    % Plot all trajectories
    figure; %hold on; grid on; axis equal;

    for j = 1:length(stateDims)
        subplot(length(stateDims), 1, j);

        hold on; grid on;

        for i = 1:length(trajectories)
            plot(t_nom, trajectories{i}(stateDims(j), :), 'b-', 'LineWidth', 0.5);
        end
    
        plot(t_nom, x_nom(stateDims(j), :), 'k--', 'LineWidth', 2);

        xlabel('time [s]');
        ylabel(['x_{', num2str(stateDims(j)), '}'])
    end

    sgtitle('Monte Carlo Rollout Trajectories');
end


function plot_input_profiles(input_profiles, u_nom, t_nom)
    % Plot all trajectories and nominal trajectory
    figure;
    
    m = size(u_nom,1);
    for i=1:m
        subplot(m,1,i)
        hold on; grid on; 
        for j = 1:length(input_profiles)
            plot(t_nom, input_profiles{j}(i, :), 'b-', 'LineWidth', 0.5);
        end
        plot(t_nom, u_nom(i, :), 'k--', 'LineWidth', 2);
        xlabel('time'); ylabel(sprintf('u_{%d}', i))

        if i==1
            title('Input history from Monte Carlo Rollouts');
        end
    end
end

% Plot error and cost metrics over time
function plot_error_metrics(errors, time)
    
    figure; hold on
     
    % Calculate upper and lower bounds of data spread
    upper_bound = mean(errors, 1) + std(errors,1);
    lower_bound = mean(errors, 1) - std(errors,1);

    % Plot the shaded std dev region
    fill([time, fliplr(time)], [upper_bound, fliplr(lower_bound)], ...
        [0.8, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Shaded region
    
    %Plot the mean
    plot(time, mean(errors, 1), 'r-', 'LineWidth', 2);
    grid on;
    xlim([0 time(end)+0.5]);
    xlabel('Time (s)'); ylabel('Error Norm');
    title('Mean of Error Norm Over Time');  
end