T = 6000;
Ns = [10 50 100 500 1000];

%N= 10; T = 30; spikes = spnet(['./sims/izhik_T_' num2str(T) '_N_10_10.mat'], T, 10);
%spnet(['./sims/izhik_T_' num2str(T) '_N_50.mat'], T, 50);
%spnet(['./sims/izhik_T_' num2str(T) '_N_100.mat'], T, 100);
%spnet(['./sims/izhik_T_' num2str(T) '_N_500.mat'], T, 500);
%spnet(['./sims/izhik_T_' num2str(T) '_N_1000.mat'], T, 1000);

parfor idx = 1:length(Ns)
	N = Ns(idx);
	spnet(['./sims/izhik_T_' num2str(T) '_' num2str(N) '_10.mat'], T, N);
end
