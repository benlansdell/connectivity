function logli = ll_network(model, data, distr)
	%Compute log likelihood of a set of data points given a set of fitted coefficients.
	%
	%
	%Usage:
	%	logli = ll(model, data)
	%     
	%Input:
	%	model = a structure of fit coefficients from MLE_glmfit
	%	data = a structure of stimulus and spike history data from ./models
	%	distr = one of 'normal', or 'poisson'
	%   
	%Output:
	%	logli = a vector of deviance for each unit
	%  
	%Test code:
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	const = 'on';
	%	nK_sp = 50; 
	%	nK_pos = 10;
	%	dt_sp = 0.002;
	%	dt_pos = 0.05;
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%	model = MLE_glmfit(data, const);
	%   logli = ll(model, data, 'poisson')

	nU = size(data.y,1);
	N = size(data.X,1);
	logli = zeros(1,nU);
	for idx = 1:nU
		%Compute y
		y = data.y(idx,:)';
		%Compute mu (using glmval)
		b_hat = model.b_hat(idx,:);
		if strcmp(distr, 'poisson')		
			mu = glmval(b_hat', data.X, 'log');
			logli(idx) = sum(y.*log(mu))-sum(mu);
		elseif strcmp(distr, 'normal')
			mu = glmval(b_hat', data.X, 'identity');
			sigma = model.sigma(idx);
			logli(idx) = -0.5*(N*log(2*pi*sigma^2)+sum(((y-mu).^2)/sigma^2));
		end
	end

