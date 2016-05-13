function [F, Q, mu] = fit_AR_LS_multi(data, R, P, d)
	%Fit an auto-regressive model to torque data of a given order using least squares.
	%
	%Usage:
	%	[F, Q, mu] = fit_AR_LS_multi(data, R, P, d)
	%     
	%Input:
	%	data = [N x 2] matrix where N is the number of data points. Could be cursor position data, 
	%			or cursor (direction, velocity) data...
	%	R = (optional, default = 1) Order of AR model to fit
	%	P = (optional, default = 1) Number of steps into the future to predict
	%	d = (optional, default = 1) Size of the steps to take to evaluate trajectories
	%   
	%Output:
	%	F = square (2P+2d(R-1)) matrix with fit coefficients
	%	Q = square covariance (2P+2d(R-1)) matrix giving errors
	%	mu = (2P+2d(R-1)) mean data
	%
	%Test code:
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	R = 5;
	%	P = 3;
	%	d = 25;
	%	data = pre.processed.torque;
	%	[F, Q, mu] = fit_AR_LS_multi(data, R, P, d);

	if (nargin < 2) R = 1; end
	if (nargin < 3) P = 1; end
	if (nargin < 4) d = 1; end
	F = zeros(2*P+2*d*(R-1));
	F((2*P+1):end,(2*P-1):(end-2)) = eye(2*d*(R-1));
	mu = zeros(2*P+2*d*(R-1),1);
	Q = zeros(2*P+2*d*(R-1));
	for p = 1:P
		residuals = [];
		%Set up matrices
		c = data(d*R:d:(end-d*(p-1)-1),1);
		r = data(d*R:-d:1,1);
		X_RU = toeplitz(c,r);
		c = data(d*R:d:(end-d*(p-1)-1),2);
		r = data(d*R:-d:1,2);
		X_FE = toeplitz(c,r);
		%Intersperse columns
		m = size(X_RU,1);
		X = reshape([X_RU;X_FE],m,[]);
		%Fit a constant term
		X = [ones(m,1),X];
		%Fit x and y parts separately (idx)
		for idx = 1:2
			Y = data((d*R+1+d*(p-1)):d:end,idx);
			%Do the least squares fit:
			beta_hat = (transpose(X)*X)\(transpose(X)*Y);
			%Save constant term
			mu(2*(p-1)+idx) = beta_hat(1);
			%save RU component
			F(2*(p-1)+idx,(2*P-2+1):2*d:end) = beta_hat(2:2:(end-1));
			%save FE component
			F(2*(p-1)+idx,(2*P-2+2):2*d:end) = beta_hat(3:2:end);
			%save for computing the covariance matrix
			residuals(:,idx) = Y-X*beta_hat;
		end
		%Compute the variance/covariance of residuals, Q.
		Q((2*(p-1)+1):(2*(p-1)+2),(2*(p-1)+1):(2*(p-1)+2)) = cov(residuals(:,1), residuals(:,2));
	end
	Q = sparse(Q);
	F = sparse(F);
	mu = sparse(mu);
	