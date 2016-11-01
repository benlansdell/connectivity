function results = run_granger(fn_in, fn_out)
	if (nargin < 2) fn_out = ''; end 

	%Load spikes.times
	load(fn_in);

	%Prepare preprocessed structure from spike times
	binsize = 0.002; %in seconds
	const = 'on';
	nK = 20; %length of spike history filter, in units of binsize
	nK_pos = 0;
	pval = 0.001;
	
	%Bin spike times
	N = length(spikes.times);
	sim_time = spikes.T;
	
	T = round(sim_time/binsize);
	
	S = zeros(N,T);
	for n = 1:N
		for idx = 1:length(spikes.times{n})
			t = max(1,round(spikes.times{n}(idx)/binsize));
			if t <= T
				S(n,t) = S(n,t) + 1;
			end
		end
	end
	
	%Load test preprocessed data
	processed.binnedspikes = S'; 
	processed.rates = processed.binnedspikes/binsize;
	processed.binsize = binsize;
	processed.tspks = spikes.times;
	processed.unitnames = {};
	processed.torque = zeros(T, 2);
	processed.dtorque = zeros(T, 2);
	processed.ddtorque = zeros(T, 2);
	for n = 1:N
	    processed.unitnames{n} = num2str(n);
	end
	data = filters_sprc_pos_network(processed, nK, nK_pos);

	fn_out = 0;
	tic
	[GCdev, GCpval, GCsig] = granger(processed, data, fn_out, pval);
	toc
	
	results.weights = GCdev;
	results.CI = zeros(size(GCdev));
	results.conn = GCsig;
	results.pvals = GCpval;

	if length(fn_out) > 0
		save(fn_out, 'results');
	end
end