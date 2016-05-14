%Load spike times and connectivity matrix from Guillaume's code
connectivity_path = './connectivity_net_sim_3_groups_rate_drive_4.mat';
spiking_path = './output_V4/'

%load spike times

%bin spike times

binsize = 0.005; %in second

%Set up S, H

lambda = 1;
alpha = 1;
gamma = 1;

%Run!
[Y, Ds] = ADMM_latent_sparse(S, H, rho, lambda, alpha, gamma)

%Compare connectivity matrix 

load(connectivity_path)