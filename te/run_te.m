function results = run_te(fn_in, fn_out)
	if (nargin < 2) fn_out = ''; end 

	%Load spikes.times
	load(fn_in);

	times = spikes.times;
	N = length(times);
	T = 1000*spikes.T;
	times{end+1} = [];
	times{end+1} = [N, T];
	for n = 1:N
		times{n} = floor(times{n}*1000);
	end

	tic
	[peakTE, CI, TEdelays_big] = ASDFTE(times, 1:30); % Now it has delay of 1ms to 30ms
	toc
	
	results.weights = peakTE;
	results.CI = CI;
	results.conn = peakTE > 1e-6;
	results.pvals = zeros(size(peakTE));

	if length(fn_out) > 0
		save(fn_out, 'results');
	end
end