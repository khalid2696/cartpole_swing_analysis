clc; clearvars; close all;

% Parameters
n = 2; % Dimensions
N = 4; % Number of ellipses
colors = lines(N);

figure; hold on; grid on; axis equal;
all_min_eigs = zeros(N, 1);

% 1. Generate and Plot random ellipses
for i = 1:N
    % Generate a random SPD matrix P
    [Q, ~] = qr(randn(n));
    L = diag(rand(n,1) * 5 + 1); 
    P = Q * L * Q'; % P = R^T * Lambda * R
    
    % Store the minimum eigenvalue
    all_min_eigs(i) = min(eig(P));
    
    % Plot the ellipse
    plot_ellipse(P, colors(i,:), sprintf('Ellipse %d', i));
end

% 2. Construct the Safe Bounding Ellipsoid (Sphere)
lambda_star = min(all_min_eigs);
P_safe = lambda_star * eye(n);

% 3. Plot the Safe Bound
plot_ellipse(P_safe, [0 0 0], 'Safe Bounding Sphere');
legend('show');
title('Outer Bound via Minimum Eigenvalue');

function plot_ellipse(P, color, name)
    theta = linspace(0, 2*pi, 100);
    circle = [cos(theta); sin(theta)];
    
    % Transform circle to ellipse: x^T P x = 1 -> x = P^(-1/2) * circle
    % We use the matrix square root
    ellipse_pts = P^(-1/2) * circle;
    
    plot(ellipse_pts(1,:), ellipse_pts(2,:), 'Color', color, ...
        'LineWidth', 2, 'DisplayName', name);
end