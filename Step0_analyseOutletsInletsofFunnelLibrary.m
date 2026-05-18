clc; clearvars; close all

load('./precomputedData/library.mat');

x_positions = [-3 -2 -1 1 2 3]; %of the funnel library
numFunnels = numel(x_positions);

%% Funnels at the bottom
level = -1; %-1 corresponds to level-down (hanging)

inlets_at_bottom = cell(numFunnels,1); outlets_at_bottom = cell(numFunnels,1);
for i = 1:length(x_positions) 
    tempDictionaryKey = [x_positions(i) level];
    funnel = funnelLibrary(num2str(tempDictionaryKey));

    inlets_at_bottom{i} = funnel.invarianceCertificates(:,:,1);
    outlets_at_bottom{i} = funnel.invarianceCertificates(:,:,end);
end

%% Funnels at the top
level = 1; %1 corresponds to level-up (upright)

inlets_at_top = cell(numFunnels,1); outlets_at_top = cell(numFunnels,1);
for i = 1:length(x_positions) 
    tempDictionaryKey = [x_positions(i) level];
    funnel = funnelLibrary(num2str(tempDictionaryKey));

    inlets_at_top{i} = funnel.invarianceCertificates(:,:,1);
    outlets_at_top{i} = funnel.invarianceCertificates(:,:,end);
end

%% Compute the largest outlets (set of all possible final states) and 
%  smallest inlets (fully contained within the inlets) -- both at top and bottom

% Relevant for swing-up maneuver
bounding_outlet_at_bottom = compute_bounding_ellipsoid(outlets_at_bottom) %contains the union of outlets (and hence each outlet) 
contained_inlet_at_top = compute_contained_ellipsoid(inlets_at_top) %is contained within intersection of inlet (and hence each inlet)

% Relevant for swing-down maneuver
bounding_outlet_at_top = compute_bounding_ellipsoid(outlets_at_top) %contains the union of outlets (and hence each outlet) 
contained_inlet_at_bottom = compute_contained_ellipsoid(inlets_at_bottom) %is contained within intersection of inlet (and hence each inlet)

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