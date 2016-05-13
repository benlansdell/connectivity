function y = AR_predict(data, F, mu, fn_out)
	%Predict AR model fit using fit_AR given a set of data
	%
	%	
	%
	%Usage:
	%	AR_predict(data, F, mu, fn_out)
	%     
	%Input:
	%	data = [N x 2] matrix where N is the number of data points. Could be cursor position data, 
	%			or cursor (direction, velocity) data...
	%	F = 
	%	mu = 
	%	fn_out = filename to save plot to
	%   
	%Test code:
	%	%Load test preprocessed data
	%	fn_out = './worksheets/diagnostics/AR_test.eps';
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	order = 1;
	%	horizon = 1;
	%	data = pre.processed.torque;
	%	[F, Q, mu] = fit_AR(data, order, horizon);
	%	AR_predict(data, F, mu, fn_out);

	clf;
	horizon = size(F,1);
	order = size(F,2)/2;
	for idx = 1:2
		subplot(2,1,idx)
		Y = data(order+1:end,idx);
		%Set up matrices
		c = data(order:end-1, 1);
		r = data(order:-1:1,1);
		X_RU = toeplitz(c,r);
		c = data(order:end-1,2);
		r = data(order:-1:1,2);
		X_FE = toeplitz(c,r);
		%Fit a constant term
		X = [ones(size(X_RU,1),1), X_RU, X_FE];
		beta_hat = [mu(idx); F(1,1:2*order,idx)];
		Y_hat = X*beta_hat;
		plot(Y);
		hold on
		plot(Y_hat, '.r');
		ylabel(['Axis ' num2str(idx)]);
		xlim([0 100])
		SSR = sum((Y-Y_hat).^2)
	end
	saveplot(gcf, fn_out);