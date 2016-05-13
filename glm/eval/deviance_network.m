function dev = deviance_network(model, data)
	%Compute deviance of a set of data points given a set of fitted coefficients. The deviance is given by:
	%
	%	D(y,mu) = 2\sum y_i ln (y_i / \mu_i) - y_i + \mu_i
	%
	%where
	%
	%	\mu_i = e^eta_i = e^{linear predictor}
	%
	%Usage:
	%	dev = deviance_network(model, data)
	%     
	%Input:
	%	model = a structure of fit coefficients from MLE_glmfit
	%	data = a structure of stimulus and spike history data from ./models
	%   
	%Output:
	%	dev = a vector of deviance for each unit
	%  
	%Test code:
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	const = 'on';
	%	nK_sp = 10; 
	%	nK_pos = 10;
	%	dt_sp = 0.002;
	%	dt_pos = 0.05;
	%	data = filters_sp_pos_network(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%	model = MLE_glmfit_network(data, const);
	%   dev = deviance_network(model, data)

	nU = size(data.y,1);
	N = size(data.y,2);
	dev = zeros(1,nU);
	for idx = 1:nU
		%Compute y
		y = data.y(idx,:)';
		%Compute mu (using glmval)
		b_hat = model.b_hat(idx,:);
		mu = glmval(b_hat', data.X(:,:), 'log');
		%Or compute mu (manually, these should be the same)
		%X = [ones(N, 1), squeeze(data.X(idx,:,:));];
		%mu = exp(X*b_hat');
		%We assume ylog(y) = 0 for y = 0, so let's ignore values where y == 0
		nz = y > 0;
		%Return deviance
		dev(idx) = 2*sum(y(nz).*log(y(nz)./mu(nz)))-2*sum(y-mu);
	end

