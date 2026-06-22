%% Script for saving it into dictionary format from raw files: saveToDictionary

clearvars; close all; clc;

load('./library.mat');

%% Load the precomputed files: swing-up/swing-down

load('./precomputedData/swing_up/nominal_trajectory_and_input.mat');

N = length(t_nom); %number of time instances (knot points)
n_x = size(x_nom,1); %dimensionality of state space
n_u = size(u_nom,1); %dimensionality of input space

funnel = struct('time_instances',NaN(1,N), 'trajectory',NaN(n_x,N), 'invarianceCertificates', NaN(n_x,n_x,N), ...
                'nominalControl', NaN(n_u,N), 'minmaxBounds', [NaN NaN]); %'numRollouts', 'controlLawFnHandle', 'controlLawParameters'

%save the time
funnel.time_instances = t_nom;

%save the nominal trajectory and input
funnel.trajectory = x_nom;
funnel.nominalControl = u_nom;

%parametric control
funnel.controlLawFnHandle = control_law_fn_handle; %API call: funnel.controlLawFnHandle(time: scalar, state: 4x1 vector)
funnel.controlLawParameters = params;

funnel.minmaxBounds = minmaxBounds;
funnel.numRollouts = numSamples;

%save the ellipsoidal certificates of invariance (NaN for a trajectory)
for k=1:N
    funnel.invarianceCertificates(:,:,k) = NaN;
end

if norm(x_nom(3,end) - pi) < 1e-3
    key = [0, 1]; %swing up
elseif norm(x_nom(3,end) - 0) < 1e-3
    key = [0, -1]; %swing-down
else
    error('Cannot handle other trajectory types');
end

funnelLibrary(num2str(key)) = funnel;

save('./library.mat','funnelLibrary');