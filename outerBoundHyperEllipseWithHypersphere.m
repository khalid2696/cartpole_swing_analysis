clc; clearvars; close all;

% Parameters
n = 2; % Dimensions
N = 4; % Number of ellipses
colors = lines(N);

figure; hold on; grid on; axis equal;

ellipsoids = cell(N,1);
% 1. Generate and plot random ellipses
for i = 1:N
    % Generate a random SPD matrix P
    [Q, ~] = qr(randn(n));
    L = diag(rand(n,1) * 5 + 1); 
    P = Q * L * Q'; % P = R^T * Lambda * R
    
    ellipsoids{i} = P;
    
    plot_ellipse(ellipsoids{i}, colors(i,:), sprintf('Ellipse %d', i));
end

% 2. Construct the Bounding Ellipsoid (Sphere) and Plot
P_bounding = compute_bounding_ellipsoid(ellipsoids);
plot_ellipse(P_bounding, [0 0 0], 'Bounding Sphere');

% 3. Construct the Contained Ellipsoid (Sphere) and Plot
P_bounding = compute_contained_ellipsoid(ellipsoids);
plot_ellipse(P_bounding, [0 1 0], 'Contained Sphere');

legend('show');
title('Outer Bound and Inner Bound via Minimum and Maximum Eigenvalues');

%% Function definitions
function P_bounding = compute_bounding_ellipsoid(ellipsoids)
    
    N = numel(ellipsoids);
    n = size(ellipsoids{1},1); %assuming all ellipsoids of the same dimensionality
    all_min_eigs = zeros(N, 1);

    for i = 1:N
        tempEllipsoidMatrix = ellipsoids{i}; 
        all_min_eigs(i) = min(eig(tempEllipsoidMatrix));
    end

    lambda_star = min(all_min_eigs);
    P_bounding = lambda_star * eye(n);
end

function P_bounding = compute_contained_ellipsoid(ellipsoids)
    
    N = numel(ellipsoids);
    n = size(ellipsoids{1},1); %assuming all ellipsoids of the same dimensionality
    all_max_eigs = zeros(N, 1);

    for i = 1:N
        tempEllipsoidMatrix = ellipsoids{i}; 
        all_max_eigs(i) = max(eig(tempEllipsoidMatrix));
    end

    lambda_star = max(all_max_eigs);
    P_bounding = lambda_star * eye(n);
end

function plot_ellipse(P, color, name)
    theta = linspace(0, 2*pi, 100);
    circle = [cos(theta); sin(theta)];
    
    % Transform circle to ellipse: x^T P x = 1 -> x = P^(-1/2) * circle
    % We use the matrix square root
    ellipse_pts = P^(-1/2) * circle;
    
    plot(ellipse_pts(1,:), ellipse_pts(2,:), 'Color', color, ...
        'LineWidth', 2, 'DisplayName', name);
end