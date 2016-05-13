function gof_timerescale_ISI(model, data, processed, alpha, fn_out)
	%Plot measures of goodness of fit of data given a fit GLM. Performs:
	%	-plots time rescaled distributions of ISI times (should be exponential)
	%	-makes QQ plot of expected vs observed ISI times
	%	-computes KS test stat for goodness of fit, and Chi-2 test
	%     
	%Input:
	%	model = data structure output by function in ./fitting (containing fitted coefficients)
	%	data = data structure output by ./models containing data used for fit
	%	alpha = significance level at which to perform test
	%	fn_out = basename filename to write plots to for each unit
	%
	%Test code:
	%	const = 'on';
	%	fn_out = './worksheets/11_30_2014/plots/timescale';
	%	nK_sp = 100; 
	%	nK_pos = 1;
	%	alpha = 0.05;
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
	%	model = MLE_glmfit(data, const);
	%	gof_timerescale_ISI(model, data, pre.processed, alpha, fn_out);

	%Rescale times
	[Lambda, T] = compensator(model, data);
	nU = size(data.y,1);

	FISIs = {};
	%For each unit
	for i = 1:nU
		%Compute interspike intervals
		ISIs = diff(Lambda{i});
		%Apply distribution function to convert what are hopefully exponential distributed values
		%to uniform values
		FISIs{i} = sort(1-exp(-ISIs));
		nSp = length(ISIs);

		clf
		hold on
		n = 1:nSp;
		plot(n/nSp, FISIs{i}, n/nSp, n/nSp);
		pm = erfinv(1-alpha/2)/sqrt(nSp);
		plot(n/nSp, n/nSp+pm, 'r--', n/nSp, n/nSp-pm, 'r--')
		xlabel('theoretical quantile')
		ylabel('quantile')
		xlim([0 1])
		ylim([0 1])
		result = 'insignificant';
		above = FISIs{i} > (n'/nSp)+pm;
		below = FISIs{i} < (n'/nSp)-pm;
		display(['Unit' num2str(i)])
		if any(above) | any(below)
			result = 'significant';
		else
			result = 'insignificant';
		end
		title(['Unit ' processed.unitnames{i} '. \alpha=' num2str(alpha) ' test result: ' result])
		saveplot(gcf, [fn_out '_unit_' processed.unitnames{i} '_test_' result '.eps']);
	end
end