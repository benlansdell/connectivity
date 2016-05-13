function [F, Q, mu] = fit_AR(data, order, horizon, every)
	%Fit an auto-regressive model to torque data of a given order using least squares. That is, fits a model of the form:
	%
	%	[x_2 y_2] ~ [mu_1, mu_2] + [x_{1}, x_{0}, y_{1}, y_{0},...]x[F_11 F_21] + [\epsilon_i^1, \epsilon_i^2]
	%	[x_3 y_3]				   [x_{2}, x_{1}, y_{2}, y_{1},...] [F_12 F_22]
	%	[...    ]				   [...                           ] [F_13 F_23]
	%											  	                [...      ]
	%
	% 	Note: the form is more complicated if a horizon greater than one is specified
	%
	%Usage:
	%	[F, Q, mu] = fit_AR(data, order, horizon, every)
	%     
	%Input:
	%	data = [N x 2] matrix where N is the number of data points. Could be cursor position data, 
	%			or cursor (direction, velocity) data...
	%	order = (optional, default = 1) Order of AR model to fit
	%	horizon = (optional, default = 1) Number of steps into the future to predict
	%	every = (optional, default = 1) Size of the steps to take to predict future trajectories
	%   
	%Output:
	%	F = [horizon, order*2, 2] matrix with fit coefficients
	%	Q = covariance matrix giving errors
	%	mu = mean data
	%
	%Test code:
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	order = 1;
	%	every = 1;
	%	data = pre.processed.torque;
	%	[F, Q, mu] = fit_AR_LS(data, order);

	if (nargin < 2) order = 1; end
	if (nargin < 3) horizon = 1; end
	if (nargin < 4) every = 1; end
	nD = size(data,2);
	F = zeros(horizon, order*nD, nD);
	mu = zeros(horizon, nD);
	residuals = zeros(size(data,1)-order, size(data,2));

	for h = 1:horizon
		for idx = 1:nD
			X = [];
			Y = data(order+1:end,idx);
			%Set up matrices
			for j = 1:nD
				c = data(order:end-h,j);
				r = data(order:-1:1,j);
				X = horzcat(X, toeplitz(c,r));
			end
			%Fit a constant term
			X = horzcat(ones(size(X,1),1), X);
			%Do the least squares fit:
			beta_hat = (transpose(X)*X)\(transpose(X)*Y);
			mu(h, idx) = beta_hat(1);
			F(h,1:nD*order,idx) = beta_hat(2:end);
			residuals(:,idx) = Y - X*beta_hat;
		end
	end
	%Compute the variance/covariance of residuals, Q.

	%Q = cov(residuals(:,1), residuals(:,2));
	Q = cov(residuals);