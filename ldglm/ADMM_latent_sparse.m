function [Y, Ds] = ADMM_latent_sparse(S, H, rho, lambda, alpha, gamma)

	%Tolerances
	eps = 1e-3;
	eps_rel = 1e-3;
	eps_abs = 1e-3;
	eps_relD = 1e-3;
	eps_absD = 1e-3;

	%Initialize
	n = size(S, 1);						%Number of neurons
	k = size(H, 1)/n;					%Number of time lags in spike history
	T = size(S, 2);						%Number of time bins
	Y = log(S+1);						%Transformed spike data			(n x T)
	Z = zeros(size(Y));					%Auxiliary variable 			(n x T)
	Lambda = zeros(size(Y));			%Lagrange multiplier			(n x T)
	Ds = zeros(n, n*k);					%Stacked connectivity matrices	(n x n*k)

	r_p = 1; 
	eps_p = 0;
	r_d = 1;
	eps_d = 0;

	%Mean centering operation (note it is self-adjoint and idempotent)
	A = eye(T) - ones(T,T)/T;											%(T x T)
	K = kron(ones(T, 1), eye(n));										%(n*T x n)
	AHinv = inv(H*A*A*H' + alpha*eye(n*k)); 							%(n*k x n*k)

	while (r_p > eps_p) & (r_d > eps_d)
		%%Update Y using ADMM
		%
		%Vectorize Y for some of this computation
		Yv = reshape(Y, [], 1);
		while gradL_rho'*inv_laplaceL_rho*gradL_rho > eps
			%Gradient in matrix form
			gradL_rho = -gradloglikelihood(S, Y) + rho*Y*A - (rho*Z + rho*Ds*H*A-Lambda)*A;
			%Vectorized gradient
			gradL_rho = reshape(gradL_rho, [], 1);
			%Compute inverse Hessian
			%Diagonal part of Hessian (n*T x n*T)
			D = spdiags(reshape(-laplaceloglikelihood(S, Y), [], 1) + rho, [0], n*T, n*T);
			Dinv = inv(D);
			d = spdiags(Dinv);
			delta = reshape(d,n,T);
			inv_laplaceL_rho = Dinv + Dinv*K*inv(T/rho*eye(n) - diag(delta*ones(T,1)))*K'*Dinv;
			Yv = Yv - inv_laplaceL_rho*gradL_rho;
		end
		%Return to matrix form
		Y = reshape(Yv, n, T);

		%Update D
		E = Ds;
		Gamma = zeros(size(E));
		r_pD = 1; eps_pD = 0; r_dD = 1; eps_dD = 0;
		while (r_pD > eps_pD) & (r_dD > eps_dD)
			for ii = 1:n
				Ds(ii,:) = ((Y(ii,:)*A - Z(ii,:) + Lambda(ii,:)/rho)*A*H' + alpha*E(ii,:) - Gamma(ii,:))*AHinv;
			end
			Ep = softmax(Ds+Gamma/alpha, gamma*T/(n*rho*alpha))
			Gamma = Gamma + alpha*(D - Ep);
			r_pD = norm(Ds-Ep, 'fro');
			r_dD = alpha*norm(E-Ep, 'fro');
			eps_pD = sqrt(n^2*k)*eps_absD + eps_relD*max(norm(Ds, 'fro'), norm(Ep, 'fro'));
			eps_dD = sqrt(n^2*k)*eps_absD + eps_relD*norm(Gamma, 'fro');
			E = Ep;
		end

		%%Update Lambda
		[U, Sigma, V] = svd((Y-Ds*H)*A+Lambda/rho);
		Zp = U*softmax(Sigma, lambda*sqrt(n*T)/rho)*V';
		Lambda = Lambda + rho*(Y*A-Zp);

		%%Determine stopping criteria
		r_p = norm((Y-Ds*H)*A-Zp, 'fro');
		r_d = rho*norm((Z-Zp)*A, 'fro');
		eps_p = sqrt(n*T)*eps_abs + eps_rel*max(norm((Y-Ds*H)*A, 'fro'), norm(Zp, 'fro'));
		eps_d = sqrt(n*T)*eps_abs + eps_rel*norm(Lambda*A, 'fro');

		%Update Z
		Z = Zp;
	end 
end

%These assume an exponential non-linearity
function G = gradloglikelihood(S, Y)
	%Luckily this can be done element-wise since each s component
	%is conditionally independent of the others given y.
	G = S - exp(Y);
end

function L = laplaceloglikelihood(S, Y)
	L = -exp(Y);
end

function y = softmax(X, thresh)
	y = sign(X).*max(0, abs(X)-thresh);
end