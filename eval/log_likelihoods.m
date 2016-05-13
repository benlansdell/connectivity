function ll = log_likelihoods(model, data, processed)
	%Compute likelihood of a fitted model given some data for each unit
	%     
	%Input:
	%	model = data structure output by function in ./fitting (containing fitted coefficients)
	%	data = data structure output by ./models containing data used for fit
	%	processed = data structure output by ./preprocess containing processed raw data
	%
	%Output:
	%	ll = log-likelihood of fitted model for each unit
	%
	%Test code:
	%	const = 'on';
	%	fn_out = './worksheets/glm_ls/20130117SpankyUtah001';
	%	nK_sp = 100; 
	%	nK_pos = 0;
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
	%	%Fit model
	%	model = MLE_glmfit(data, const);
	%	%Plot predictions
	%	ll = log_likelihoods(model, data, pre.processed);

	nU = size(data.y,1); %number of units
	ll = zeros(nU,1);

	%For each unit, predicting firing rate, smooth and compare to actual smoothed firing rate
	for idx=1:nU 
		b_hat = model.b_hat(idx,:);
		% Uses the fit coefficients and the original input data to generate the ouput rho
		rho_hat = glmval(b_hat', squeeze(data.X(idx,:,:)), 'log');
		%Compute likelihood
		l = log_likelihood(rho_hat', data.y(idx,:));
		ll(idx,1) = l;
	end
end