T = 6000;
Ns = [10 50 100 500 1000];
ns = [10, 20, 50, 100];

for N = Ns
	fn_in = ['./sims/izhik_T_' num2str(T) '_' num2str(N) '_10.mat'];
	for n = ns
		if n <= N
			fn_out = ['./results/te_izhik_T_' num2str(T) '_' num2str(N) '_10_n_' num2str(n) '.mat'];
			run_te(fn_in, fn_out, n);
		end
	end
end