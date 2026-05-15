load('../precomputedData/nominal_trajectory_and_input.mat');
load('../precomputedData/TVLQR_gains_and_cost_matrices.mat');

N = 150; %specify the resolution of output

N_fine = numel(t_nom);

indices = linspace(1,N_fine,N);

n_x = size(x_nom,1);
n_u = size(u_nom,1);

t_coarse = NaN(1, N);
x_coarse = NaN(n_x,N);
u_coarse = NaN(n_u,N);
K_coarse = NaN(n_u,n_x,N);
P_coarse = NaN(n_x,n_x,N);

for k = 1:N
    tempIndex = round(indices(k));
    t_coarse(k) = t_nom(tempIndex);
    x_coarse(:,k) = x_nom(:,tempIndex);
    u_coarse(:,k) = u_nom(:,tempIndex);
    K_coarse(:,:,k) = tvlqr.K(:,:,tempIndex);
    P_coarse(:,:,k) = tvlqr.S(:,:,tempIndex);
end

%% Some debug print statements to verify the downsampling
% disp(t_coarse(1)); disp(t_nom(1)); fprintf('\n'); 
% disp(x_coarse(:,end)); disp(x_nom(:,end)); fprintf('\n');  
% disp(u_coarse(:,end)); disp(u_nom(:,end)); fprintf('\n'); 
% disp(t_coarse(end)); disp(t_nom(end)); fprintf('\n'); 
% 
% disp(P_coarse(:,:,1)); disp(tvlqr.S(:,:,1)); fprintf('\n'); 
% disp(K_coarse(:,:,end)); disp(tvlqr.K(:,:,end)); fprintf('\n');

%% Saving up the data (using computeSOSFunnel repo naming convention)

time_instances = t_coarse;
x_nom = x_coarse;
u_nom = u_coarse;
K = K_coarse;
P = P_coarse;

save('./nominalTrajectory.mat',"time_instances","x_nom","u_nom");
save('./LQRGainsAndCostMatrices.mat',"time_instances","K","P");
