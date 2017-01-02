function spikes = spnet_ld_input(fn_out, T, N, n)
	% spnet.m: Spiking network with axonal conduction delays and STDP
	% Created by Eugene M.Izhikevich.                February 3, 2004
	% Modified to allow arbitrary delay distributions.  April 16,2008

	%Modified to turn off STDP halfway through....

	if (nargin < 1) fn_out = ''; end
	if (nargin < 2) T = 2; end 					%Test duration
	if (nargin < 3) N = 1000; end 
	if (nargin < 4) n = 10; end 
	if (N < 10) N = 10; end 

	display(['Simulating ' num2str(N) ' neurons'])

	rand('seed',1);
	M=max(ceil(0.1*N),2);  % number of synapses per neuron
	D=20;                  % maximal conduction delay 
	% excitatory neurons   % inhibitory neurons      % total number 
	Ne=ceil(0.8*N);                Ni=N-Ne; 
	a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
	d=[   8*ones(Ne,1);    2*ones(Ni,1)];
	sm=20;                 % maximal synaptic strength

	exw = 6;
	inw = 5;

	%Low-dim input details
	nC = 3;
	L = 6;
	ep = 1;
	ep_C = 5;
	%%Choose some random stable dynamical system
	A = rand(L,L,nC);
	for idx = 1:nC
		[V,lambda] = eig(A(:,:,idx));
		lambda = diag(lambda);
		A(:,:,idx) = A(:,:,idx)/(max(lambda)+.2);
	end

	%Initialize 
	x = zeros(L,1);

	%Generate loading vectors
	%Choose n indices that won't be stimulated
	recorded = randsample(N,n);
	C = ep_C*rand(N,L,nC);
	C(recorded,:,:) = 0;

	%Generate connectivity and weights matrix
	conn = zeros(N,N);
	weight = zeros(N,N);
	delaysmat = zeros(N,N);

	% post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
	% Take special care not to have multiple connections between neurons
	delays = cell(N,D);
	for i=1:Ne
		p=randperm(N);
		post(i,:)=p(1:M);
		for j=1:M
			delay = ceil(D*rand);
			delays{i, delay}(end+1) = j;  % Assign random exc delays
		end
	end

	for i=Ne+1:N
		p=randperm(Ne);
		post(i,:)=p(1:M);
		delays{i,1}=1:M;                    % all inh delays are 1 ms.
	end
	
	s=[exw*ones(Ne,M);-inw*ones(Ni,M)];         % synaptic weights
	sd=zeros(N,M);                          % their derivatives
	
	%Only needed for STDP
	% Make links at postsynaptic targets to the presynaptic weights
	pre = cell(N,1);
	aux = cell(N,1);
	for i=1:Ne
		for j=1:D
			for k=1:length(delays{i,j})
				pre{post(i, delays{i, j}(k))}(end+1) = N*(delays{i, j}(k)-1)+i;
				aux{post(i, delays{i, j}(k))}(end+1) = N*(D-1-j)+i; % takes into account delay
			end
		end
	end
	
	STDP = zeros(N,1001+D);
	v = -65*ones(N,1);                      % initial values
	u = 0.2.*v;                             % initial values
	firings=[-D 0];                         % spike timings
	
	seco = 0;
	tic;
	display(['Simulating for ' num2str(T) ' seconds with STDP'])
	for sec=1:T                      		% simulation of T seconds
		regime = ceil(nC*sec/T);
		%Progress 
		if (floor(sec/T*10) ~= floor(seco/T*10)) display([num2str(floor(sec/T*100)) '%']); end;
		seco = sec;
		for t=1:1000                          % simulation of 1 sec

			%I=zeros(N,1);
			%I(ceil(N*rand))=20;                 % random thalamic input 

			%Update LD dynamics 
			x = A(:,:,regime)*x + ep*randn(L,1);
			%Generate thalamic input 
			I = max(C(:,:,regime)*x,0);

			fired = find(v>=30);                % indices of fired neurons
			v(fired)=-65;  
			u(fired)=u(fired)+d(fired);
			STDP(fired,t+D)=0.1;
			for k=1:length(fired)
				sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
			end;
			firings=[firings;t*ones(length(fired),1),fired];
			k=size(firings,1);
			while firings(k,1)>t-D
				del=delays{firings(k,2),t-firings(k,1)+1};
				ind = post(firings(k,2),del);
				I(ind)=I(ind)+s(firings(k,2), del)';
				sd(firings(k,2),del)=sd(firings(k,2),del)-1.2*STDP(ind,t+D)';
				k=k-1;
			end
			v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical 
			v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time 
			u=u+a.*(0.2*v-u);                   % step is 0.5 ms
			STDP(:,t+D+1)=0.95*STDP(:,t+D);     % tau = 20 ms
		end
		%plot(firings(:,1),firings(:,2),'.');
		%axis([0 1000 0 N]); drawnow;
		STDP(:,1:D+1)=STDP(:,1001:1001+D);
		ind = find(firings(:,1) > 1001-D);
		firings=[-D 0; firings(ind,1)-1000,firings(ind,2)];

		%Comment out to ignore STDP
		s(1:Ne,:)=max(0,min(sm,0.01+s(1:Ne,:)+sd(1:Ne,:)));
		sd=0.9*sd;
	end

	%Turn off STDP and note weights and connections
	for i=1:Ne
		conn(i,post(i,:)) = 1;
		weight(i,post(i,:)) = s(i,:);
		for del=1:D
			units = delays{i,del};
			if length(units) > 0
				delaysmat(i,post(i,units)) = del;
			end
		end;
	end;
	for i=Ne+1:N
		conn(i,post(i,:)) = 1;
		weight(i,post(i,:)) = s(i,:);
		delaysmat(i,post(i,:)) = 1;
	end;

	times = cell(N,1);
	seco = 0;
	display(['Simulating for ' num2str(T) ' seconds without STDP'])
	for sec=1:T                      		% simulation of T seconds
		regime = ceil(nC*sec/T);
		%Progress 
		if (floor(sec/T*10) ~= floor(seco/T*10)) display([num2str(floor(sec/T*100)) '%']); end;
		seco = sec;
		for t=1:1000                          % simulation of 1 sec
			%I=zeros(N,1);
			%I(ceil(N*rand))=20;                 % random thalamic input 

			%Update LD dynamics 
			x = A(:,:,regime)*x + ep*randn(L,1);
			%Generate thalamic input 
			I = max(C(:,:,regime)*x,0);

			fired = find(v>=30);                % indices of fired neurons
			if length(fired) > 0 
				for i = fired'
					times{i}(end+1) = (sec-1)+t/1000;
				end
			end
			v(fired)=-65;  
			u(fired)=u(fired)+d(fired);
			firings=[firings;t*ones(length(fired),1),fired];
			k=size(firings,1);
			while firings(k,1)>t-D
				del=delays{firings(k,2),t-firings(k,1)+1};
				ind = post(firings(k,2),del);
				I(ind)=I(ind)+s(firings(k,2), del)';
				k=k-1;
			end;
			v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical 
			v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time 
			u=u+a.*(0.2*v-u);                   % step is 0.5 ms
		end;
		%plot(firings(:,1),firings(:,2),'.');
		%axis([0 1000 0 N]); drawnow;
		ind = find(firings(:,1) > 1001-D);
		firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
	end
	t = toc;
	display(['Time elapsed ' num2str(t) ' seconds.'])

	spikes.times = times;		%Spike times for each unit
	spikes.units = '1'; 		%Units in seconds
	spikes.conn = conn;			%Connectivity matrix
	spikes.weights = weight;	%Connection weight (negative for inhibitory)
	spikes.delays = delaysmat;	%Delay matrix (if applicable)
	spikes.N = N;				%Number of units
	spikes.T = T;
	spikes.A = A;				%Dynamical systems (transition matrices)
	spikes.C = C;				%Loading matrix (stays constant for different regimes)
	%spikes.input = inp;		%Input vector

	%Save results
	if length(fn_out) > 0
		save(fn_out, 'spikes');
	end
end