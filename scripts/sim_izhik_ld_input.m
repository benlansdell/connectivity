T = 6000;
Ns = [10 50 100];
ns = [5 10 15];

parfor idx = 1:length(Ns)
	for j = 1:length(ns)
		n = ns(j);
		N = Ns(idx);
		fn_out = ['./sims/ld_input_izhik_T_' num2str(T) '_' num2str(N) '_10.mat'];
		if (n < N) && exist(fn_out, 'file') ~= 2
			spnet_ld_input(fn_out, T, N, n);
		end
	end
end