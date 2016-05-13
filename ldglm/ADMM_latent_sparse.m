function [Y, D] = ADMM_latent_sparse(S, H, rho, lambda, alpha, gamma)

	%Tolerances
	eps = 1e-3;
	eps_rel = 1e-3;
	eps_abs = 1e-3;
	eps_relD = 1e-3;
	eps_absD = 1e-3;

	%Initialize
	n = size(S, 1);
	T = size(S, 2);
	Y = log(S+1);
	Z = zeros(size(Y));
	Lambda = zeros(size(Y));
	D = zeros(n);
	r_p = 1; 
	eps_p = 0;
	r_d = 1;
	eps_d = 0;

	%Mean centering operation (note it is self-adjoint and idempotent)
	A = eye(T) - ones(T,T)/T;
	K = kron(ones(T, 1), eye(n));
	AHinv = inv(H*A*A*H' + alpha*eye(n*k));

	while (r_p > eps_p) & (r_d > eps_d)
		%%Update Y
		%
		%Vectorize Y for some of this computation
		Yv = reshape(Y, [], 1);
		while gradL_rho'*inv_laplaceL_rho*gradL_rho > eps
			%Gradient in matrix form
			gradL_rho = -gradloglikelihood(S, Y) + rho*Y*A - (rho*Z + rho*D*H*A-Lambda)*A;
			%Vectorized gradient
			gradL_rho = reshape(gradL_rho, [], 1);
			%Diagonal part of Hessian (nT x nT)
			D = spdiags(reshape(-laplaceloglikelihood(S, Y), [], 1) + rho);
			Dinv = inv(D);
			d = spdiags(Dinv);
			delta = reshape(d,n,T);
			inv_laplaceL_rho = Dinv + Dinv*K*(T/rho*eye(n) - diag(delta*ones(T,1)))*K'*Dinv;
			Yv = Yv - inv_laplaceL_rho*gradL_rho;
		end
		%Return to matrix form
		Y = reshape(Yv, n, T);

		%Update D
		E = D;
		Gamma = zeros(size(E));
		r_pD = 1; eps_pD = 0; r_dD = 1; eps_dD = 0;
		while (r_pD > eps_pD) & (r_dD > eps_dD)
			for ii = 1:n
				D(ii,:) = ((Y(ii,:)*A - Z(ii,:) + Lambda(ii,:)/rho)*A*H' + alpha*E(ii,:) - Gamma(ii,:))*AHinv;
			end
			Ep = softmax(D+Gamma/alpha, gamma*T/(n*rho*alpha))
			Gamma = Gamma + alpha*(D - Ep);
			r_pD = norm(D-Ep, 'fro');
			r_dD = alpha*norm(E-Ep, 'fro');
			eps_pD = sqrt(n^2*k)*eps_absD + eps_relD*max(norm(D, 'fro'), norm(Ep, 'fro'));
			eps_dD = sqrt(n^2*k)*eps_absD + eps_relD*norm(Gamma, 'fro');
			E = Ep;
		end

		%%Update Lambda
		[U, Sigma, V] = svd((Y-D*H)*A+Lambda/rho);
		Zp = U*softmax(Sigma, lambda*sqrt(n*T)/rho)*V';
		Lambda = Lambda + rho*(Y*A-Zp);

		%%Determine stopping criteria
		r_p = norm((Y-D*H)*A-Zp, 'fro');
		r_d = rho*norm((Z-Zp)*A, 'fro');
		eps_p = sqrt(n*T)*eps_abs + eps_rel*max(norm((Y-D*H)*A, 'fro'), norm(Zp, 'fro'));
		eps_d = sqrt(n*T)*eps_abs + eps_rel*norm(Lambda*A, 'fro');

		%Update Z
		Z = Zp;
	end 
end

%These assume an exponential non-linearity
function gradloglikelihood(S, Y)
	%Luckily this can be done element-wise since each s component
	%is conditionally independent of the others given y.
	G = S - exp(Y);
end

function laplaceloglikelihood(S, Y)
	L = -exp(Y);
end

function y = softmax(X, thresh)
	y = sign(X).*max(0, abs(X)-thresh);
end