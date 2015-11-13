function [F, Q, mu] = fit_AR(data, order, horizon, every)
	%Fit an auto-regressive model to torque data of a given order using arfit package.
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
	%	[F, Q, mu] = fit_AR(data, order);

	if (nargin < 2) order = 1; end
	if (nargin < 3) horizon = 1; end
	if (nargin < 4) every = 1; end
	F = zeros(horizon, order*2, 2);
	mu = zeros(horizon, 2);
	residuals = zeros(size(data,1)-order, size(data,2));
	pmin = order;
	pmax = order;

	for h = 1:horizon
		[w, A, C] = arfit(data, pmin, pmax);
	end
	%Compute the variance/covariance of residuals, Q.
