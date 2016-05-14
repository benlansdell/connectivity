%Load spike times and connectivity matrix from Guillaume's code
connectivity_path = './connectivity_net_sim_3_groups_rate_drive_4.mat';
spiking_path = './output_V4/Net_sim_3_groups_rate_driven_4_noSTDP_spiketimes_batch_1_10sec.mat';
otherdata_path = './output_V4/Net_sim_3_groups_rate_driven_4_noSTDP_10sec.mat';

load(spiking_path);
load(otherdata_path);
load(connectivity_path);

binsize = 0.005; %in seconds
K = 6; %length of spike history filter, in units of binsize

%Bin spike times
N = length(spiketimes);
sim_time = time(2);
T = round(sim_time/binsize);

%S -- spike trains -- N x T
S = zeros(N,T);
for n = 1:N
	for idx = 1:length(spiketimes{n})
		t = max(1,round(spiketimes{n}(idx)/binsize));
		S(n,t) = S(n,t) + 1;
	end
end

%H -- spike histories -- (N*K) x T
H = zeros(N*K,T);
for k = 1:K
	H((k-1)*N+1:k*N,k:end) = S(:,1:(end-k+1));
end

rho = 1;
alpha = 1;
gamma = 1;

%Supposedly the only free parameter...
lambda = 1;

%Run!
[Y, Ds] = ADMM_latent_sparse(S, H, rho, lambda, alpha, gamma);

%Compare connectivity matrices
