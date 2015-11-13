function processed = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU, k_FE, N, binsize, seed)
	% Function to generate data according to a provided GLM. That is, given a set of filters, will simulate 
	%	filtered Gaussian white noise torque data, run that through a GLM, and produce an output spike train.
	%
	%Usage:
	%	processed = generate_glm_data(freqlow, freqhigh, k_sp, k_RU, k_FE, N, binsize)
	%
	%Input:
	%	freqlow = amplitude of low frequency input stimulus data (in Hertz)
	%	freqhigh = amplitude of high frequency (noise) input stimulus data (in Hertz)
	%	k_const = constant firing rate
	%	k_sp = spike history filter
	%	k_RU = RU cursor data filter
	%	k_FE = FE cursor data filter
	%	N = (optional, default = 10000) number of data points
	%	binsize = (optional, default = 0.002) size of bins in seconds
	%	seed = (optional) if provided will run rng(seed) to that torque data is reproducible
	%		
	%Output:
	%	processed is a structure containing the following fields:
	%		binnedspikes = [nB x nU] array with spikes from all channels binned according to binsize. nB = no. bins, nU = no. units.
	%		rates = [nB x nU] array of binnedspikes multiplied by samplerate to give a per second firing rate approximation
	%		torque = [nB x 2] array of torque inputs 
	%		dtorque = [nB x 2] array of diff of torque inputs (approximation of velocity)
	%		ddtorque = [nB x 2] array of diff of diff of torque inputs (approximation of acceleration)
	%		unitnames = String of the format Electrode.SortCode used to distinguish individual units from multi-unit electrode activity
	%		tspks = cell array containing spike times for each active channel
	%		binsize = binsize used
	%
	%Test code:
	%		freqlow = 7;
	%		freqhigh = 4000;
	%		k_const = 0;
	%		k_sp = 0;
	%		k_RU = 0;
	%		k_FE = 0;
	%		N = 10000;
	%		binsize = 0.002;
	%		processed = generate_glm_data(freqlow, freqhigh, k_sp, k_RU, k_FE, N, binsize);
	
	%Optional arguments
	if (nargin < 7) N = 10000; end
	if (nargin < 8) binsize = 0.002; end
	if (nargin == 9) rng(seed); end

	nU = 1;
	nK_sp = length(k_sp);
	nK_pos = length(k_RU);
	dt_sp = binsize;
	dt_pos = 0.05; %50ms timebins

	%Simulate then filter an OU (mean reverting) process
	%Generate random data
	kap = 1;
	gam = 1;
	rRU = randn(N, 1);
	rFE = randn(N, 1);	
	tRU = zeros(N, 1);
	tFE = zeros(N, 1);
	for idx = 2:N
		tRU(idx) = tRU(idx-1)+binsize*(-kap*tRU(idx-1)+rRU(idx))/gam;
		tFE(idx) = tFE(idx-1)+binsize*(-kap*tFE(idx-1)+rFE(idx))/gam;
	end

	%Create a stop filter that blocks any frequnecies between low and high
	n = 8;
	nyquist = 1/binsize/2;
	Wn = [freqlow, freqhigh]/nyquist;
	if (freqhigh < nyquist)
		[b, a] = butter(n,Wn,'stop');	
	else
		[b, a] = butter(n,Wn(1),'low');
	end
	torque = [filter(b, a, tRU), filter(b, a, tFE)];

	%Plot
	%plot(cumsum(tRU)); hold on
	plot(torque(:,1), torque(:,2), 'r');

	%Set up input structure
	processed.binnedspikes = zeros(N, nU);
	processed.rates = processed.binnedspikes/binsize;
	processed.torque = torque;
	processed.dtorque = diff(torque);
	processed.ddtorque = diff(processed.dtorque);
	processed.unitnames = {'gen. GLM train'};
	processed.tspks = 0;
	processed.binsize = binsize;
	%Produce data ready to apply filter to
	data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%Setup model structure
	model.b_hat = [k_const, k_sp, k_RU, k_FE];
	%Simulate GLM
	[y, tspks] = glmsim(processed, model, data);
	%Return spike train
	processed.binnedspikes = y;
	processed.rates = y/binsize;
	processed.tspks = tspks;
	processed.torque = data.torque;
	processed.dtorque = data.dtorque;
	processed.ddtorque = data.ddtorque;