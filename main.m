clc; clearvars; close all;

%parpool; %initialise parallel processing %if clusters available and toolboox installed

%% Add directories and create if they don't exist
addpath('./lib/');

%% Specify Cart-Pole Parameters (not defining these will result in an error)

cartPoleParameters.M = 1.0;    % cart mass (kg)
cartPoleParameters.m = 0.1;    % pendulum mass (kg)
cartPoleParameters.L = 0.5;    % pendulum length (m)
cartPoleParameters.g = 9.81;   % gravity (m/s^2)

%% Define the system dynamics as a function handle
dynamicsFnHandle = @(x, u) cartpole_dynamics(x, u, cartPoleParameters);

%% Compute a nominal trajectory and corresponding feedforward inputs (using swing-up control law)

run('./swingUpEnergyControlLaw_IVP.m');

keyboard
%% for debugging purposes
run('./utils/checkOpenLoop_IVP.m');
% run('./utils/checkFeedbackController.m');

% keyboard;

%% Function definitions

% Define the system dynamics: cartpole
% x = [p_x; v_x; theta; theta_dot]
% u = F
function f = cartpole_dynamics(x, u, cartPoleParameters)
    % Numerical evaluation of cartpole dynamics
    
    
    %extract parameters
    M = cartPoleParameters.M; m = cartPoleParameters.m;
    L = cartPoleParameters.L; g = cartPoleParameters.g;

    % Extract states
    p_x = x(1);
    v_x = x(2);
    theta = x(3);
    omega = x(4);
    
    % Control input
    F = u;
    
    % Define trigonometric functions
    s_theta = sin(theta);
    c_theta = cos(theta);
    
    % Common denominator
    %denom = M + m*(1 - c_theta^2);
    denom = M + m*s_theta^2;
    
    % State derivatives
    f = [
        v_x;
        (F + m*L*omega^2*s_theta + m*g*s_theta*c_theta) / denom;
        omega;
        (-F*c_theta - m*L*omega^2*s_theta*c_theta - (M + m)*g*s_theta) / (L * denom)
    ];
end
