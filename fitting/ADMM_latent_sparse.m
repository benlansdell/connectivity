function [Y, D] = ADMM_latent_sparse(S, H, rho, lambda, gamma)

	%Tolerances
	eps = 1e-3;
	eps_rel = 1e-3;
	eps_abs = 1e-3;

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

	while (r_p > eps_p) & (r_d > eps_d)
		while gradL_rho'*inv_laplaceL_rho*gradL_rho > eps
			gradL_rho = -gradloglikelihood(S, Y) + rho*rowmean(Y) - rowmean(rho*Z + rho*rowmean(D*H)-Lambda)';
			K = kron(ones(T, 1), eye(n));
			D = -laplaceloglikelihood(S, Y) + rho*eye(n*T);
			Dinv = inv(D);
			d = diag(Dinv);
			delta = reshape(d, n, T);
			inv_laplaceL_rho = Dinv + Dinv*K*(T/rho*eye(n) - diag(delta*ones(T,1)))*K'*Dinv;
			Y = Y - inv_laplaceL_rho*gradL_rho;
		end
		[U, Sigma, V] = svd(rowmean(Y)+Gamma/rho);
		Zp = U*Softmax(Sigma, lambda*sqrt(n*T)/rho)*V;
		Lambda = Lambda + rho*(rowmean(Y)-Zp);
		r_p = norm(rowmean(Y)-Zp);
		r_d = rho*norm(rowmean(Z-Zp)');
		eps_p = sqrt(n*T)*eps_abs + eps_rel*max(norm(rowmean(Y)), norm(Zp));
		eps_d = sqrt(n*T)*eps_abs + eps_rel*norm(rowmean(Lambda)');
		Z = Zp;
	end 
end

function rowmean(Y)


end

function gradloglikelihood(S, Y)
	%Luckily this can be done element-wise since each s component
	%is conditionally independent of the others given y.

end

function laplaceloglikelihood(S, Y)

end

function softmax(X, thresh)

end