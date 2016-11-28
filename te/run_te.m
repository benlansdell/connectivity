function results = run_te(fn_in, fn_out, n)
	if (nargin < 2) fn_out = ''; end 
	if (nargin < 3) n = 0; end 

	%Load spikes.times
	load(fn_in);

	times = spikes.times;
	N = length(times);
	if (n == 0) 
		J = 1:N;
		n = N;
	else
		J = randsample(1:N,n);
	end
	J = sort(J);
	T = 1000*spikes.T;

	display(['Estimating TE for ' num2str(spikes.T) ' seconds'])
	times = times(J);
	times{end+1} = [];
	times{end+1} = [n, T];
	for idx = 1:n
		times{idx} = floor(times{idx}*1000);
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