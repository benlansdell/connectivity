function processed = generate_glm_data_torque(pre, k_const, k_sp, k_RU, k_FE, dt_sp, dt_pos, N, binsize, seed)
	% Function to generate data according to a provided GLM. That is, given a set of filters, will use recorded
	% torque data, run that through a given GLM, and produce a simulated spike train.
	%
	%Usage:
	%	processed = generate_glm_data_torque(k_sp, k_RU, k_FE, N, binsize, seed)
	%
	%Input:
	%	pre = processed structure
	%	k_const = constant firing rate
	%	k_sp = spike history filter
	%	k_RU = RU cursor data filter
	%	k_FE = FE cursor data filter
	%	N = (optional, default = 10000) number of data points
	%	binsize = (optional, default = 0.002) size of bins in seconds
	%	seed = (optional) if provided will run rng(seed) so that torque data is reproducible
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
	%Example code:
	%		%Load torque data from recording
	%		nevfile = './testdata/20130117SpankyUtah001.nev';
	%		binsize = 0.002;
	%		offset = 0.0;
	%		threshold = 5;
	%		processed = preprocess(nevfile, binsize, threshold, offset, fn_out);
	%		k_const = 0;
	%		k_sp = 0;
	%		k_RU = 0;
	%		k_FE = 0;
	%		N = 10000;
	%		binsize = 0.002;
	%		%Return structure with actual spikes replaced by simulated spikes given GLM
	%		processed = generate_glm_data(processed, k_sp, k_RU, k_FE, N, binsize);
	
	%Optional arguments
	if (nargin < 6) N = 10000; end
	if (nargin < 7) binsize = 0.002; end
	if (nargin == 8) rng(seed); end

	nU = 1;
	nK_sp = length(k_sp);
	nK_pos = length(k_RU);

	%Set up input structure
	processed.binnedspikes = zeros(N, nU);
	processed.rates = processed.binnedspikes/binsize;
	processed.torque = pre.torque;
	processed.dtorque = pre.dtorque;
	processed.ddtorque = pre.ddtorque;
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
	%These have now been truncated, so are different than what's defined above
	processed.torque = data.torque;
	processed.dtorque = data.dtorque;
	processed.ddtorque = data.ddtorque;