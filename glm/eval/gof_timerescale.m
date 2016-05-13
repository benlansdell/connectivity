function gof_timerescale(model, data, processed, alpha, fn_out)
	%Plot measures of goodness of fit of data given a fit GLM. Performs:
	%	-plots time rescaled distributions of ISI times (should be exponential)
	%	-makes QQ plot of expected vs observed ISI times
	%	-computes KS test stat for goodness of fit, and Chi-2 test
	%	-deviance of model, compared to deviance expected from GLM (should be similar, if model is not overfit)
	%     
	%Input:
	%	model = data structure output by function in ./fitting (containing fitted coefficients)
	%	data = data structure output by ./models containing data used for fit
	%	alpha = significance level at which to perform test
	%	fn_out = basename filename to write plots to for each unit
	%
	%Test code:
	%	const = 'on';
	%	fn_out = './worksheets/11_16_2014/plots/timescale';
	%	nK_sp = 100; 
	%	nK_pos = 1;
	%	alpha = 0.05;
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
	%	model = MLE_glmfit(data, const);
	%	gof_timerescale(model, data, pre.processed, alpha, fn_out);

	%Rescale times
	[Lambda, T] = compensator(model, data);
	nU = size(data.y,1);

	%For each unit plot
	for i = 1:nU
		clf
		Ti = T(i);
		Ni = length(Lambda{i});
		n = 1:Ni;
		hold on
		plot(Lambda{i}/Ti, n/Ni, n/Ni, n/Ni);
		pm = erfinv(1-alpha/2)/sqrt(Ti);
		plot(n/Ni, n/Ni+pm, 'r--', n/Ni, n/Ni-pm, 'r--')
		xlabel('\tau(t)/\tau(T)')
		ylabel('N_t/N_T')
		xlim([0 1])
		ylim([0 1])
		above = n' > (Lambda{i}/Ti+pm)*Ni;
		below = n' < (Lambda{i}/Ti-pm)*Ni;
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