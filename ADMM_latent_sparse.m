function [Y, D] = ADMM_latent_sparse(S, H, rho, lambda, gamma)

	nU = size(S, 1);
	nB = size(S, 2);
	Y = log(S+1);
	Z = zeros(size(Y));
	D = zeros(nU);
	r_p = 1; 
	eps_p = 0;
	r_d = 1;
	eps_d = 0;
	eps = 1e-3;
	eps_rel = 1e-3;
	eps_abs = 1e-3;

	while (r_p > eps_p) & (r_d > eps_d)
		while gradL_rho'*laplaceL_rho*gradL_rho > eps



		end


	end 
end

function rowmean(Y)


end

function loglikelihood(S, Y)


end

