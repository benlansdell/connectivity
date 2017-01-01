%Load spiking data
T = 6000;
N = 50;
fn_in = ['./sims/izhik_T_' num2str(T) '_' num2str(N) '_10.mat'];
load(fn_in);

load(spiking_path);
load(otherdata_path);
load(connectivity_path);

binsize = 0.005; %in seconds
K = 3; %length of spike history filter, in units of binsize
ntrunc = 40;

%Truncate to a smaller set of neurons
%spiketimes = spiketimes(1:ntrunc);

%Bin spike times
N = length(spikes.times);
sim_time = spikes.T;
J = spikes.weights;

%Chop in quarter...
sim_time = sim_time/4;

T = round(sim_time/binsize);

%S -- spike trains -- N x T
S = zeros(N,T);
for n = 1:N
	for idx = 1:length(spikes.times{n})
		t = max(1,round(spikes.times{n}(idx)/binsize));
		if t <= T
			S(n,t) = S(n,t) + 1;
		end
	end
end

%H -- spike histories -- (N*K) x T
H = zeros(N*K,T);
for k = 1:K
	H((k-1)*N+1:k*N,(k+1):end) = S(:,1:(end-k));
end

rho =   1;
alpha = 1;

exps = -7:2:3;

for l = exps 
	for g = exps
		lambda = l;  	%low rank penalty
		gamma  = g;		%sparsity penalty
		[Y, Ds] = ADMM_LDGLM(S, H, rho, lambda, alpha, gamma);
		save(['./worksheets/gamma_lambda_sweep/gamma_' num2str(g) '_lambda_' num2str(l) '.mat'], 'Y', 'Ds', 'K', 'N', 'J')

		%Compare connectivity matrices
		Dtot = zeros(ntrunc, ntrunc)
		for k = 1:K
			D = Ds(:,((k-1)*N+1):(k*N));
			for ii = 1:N
				D(ii,ii) = 0;
			end
			Dtot = Dtot + abs(D);
		end
		a = reshape(J(1:N,1:N), 1, []);
		b = reshape(Dtot, 1, []);
		c = corr(a',b')^2
		figure 
		plot(a, b, '.k')
		xlabel('actual')
		ylabel('estimated')
		title(['Real vs predicted connectivity r^2: ' num2str(c)])
		saveplot(gcf, ['./worksheets/gamma_lambda_sweep/gamma_' num2str(g) '_lambda_' num2str(l) '_connectivity_.eps'])		
		%Plot Y vs spike raster for sample unit 
		clf
		sptimes = find(S(ntrunc,:)>0);
		hold on 
		for i = 1:length(sptimes)
			plot(sptimes(i), 0, 'k.')
		end
		plot(exp(Y(n,:)))
		xlabel('time')
		ylabel('rate')
		xlim([0 400])
		saveplot(gcf, ['./worksheets/gamma_lambda_sweep/gamma_' num2str(g) '_lambda_' num2str(l) '_y.eps'])
	end
end
