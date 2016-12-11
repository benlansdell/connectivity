%Load spiking data
T = 6000;
N = 50;
fn_in = ['./sims/izhik_T_' num2str(T) '_' num2str(N) '_10.mat'];
load(fn_in);

binsize = 0.05; %in seconds
K = 0; %length of spike history filter, in units of binsize
%ntrunc = 60;

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

%Supposedly the only free parameters...
lambda = 0.03;  	%low rank penalty
gamma = 1e-6;	%sparsity penalty

%But what should these be (see nucnrmmin.m for guidance)?
rho =   1.3;
alpha = 1;

%Run!
[Y, Ds] = ADMM_latent_sparse(S, H, rho, lambda, alpha, gamma);

%Also fun Pfau's version
opts.lambda = lambda;
opts.nlin = 'exp';
opts.center = true;
opts.q = k;
opts.gamma = gamma;
[Ynn, Dnn] = nucnrmmin(S, opts);

%Compare connectivity matrices

%for k = 1:K
figure 
for k = 1:K
	D = Ds(:,((k-1)*N+1):(k*N));
	for ii = 1:N
		D(ii,ii) = 0;
	end
	a = reshape(J(1:N,1:N), 1, []);
	b = reshape(D, 1, []);
	c = corr(a',b')^2
	subplot(ceil(sqrt(K)), ceil(sqrt(K)), k);
	plot(a, b, '.k')
	xlabel('actual')
	ylabel('estimated')
	title(['Real vs predicted connectivity r^2: ' num2str(c)])
end
saveplot(gcf, './connectivity.eps')

%Plot Y vs spike raster for sample unit 
figure 
clf
sptimes = find(S(:,:)>0);
hold on 
for i = 1:length(sptimes)
	plot(sptimes(i), 0, 'k.')
end
plot(exp(Y(n,:)))
xlabel('time')
ylabel('rate')
xlim([0 400])
saveplot(gcf, './sample_y.eps')
