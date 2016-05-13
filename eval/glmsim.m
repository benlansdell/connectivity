function [y, tspks, rho, dev] = glmsim(processed, model, data, maxspks)
	%Compute deviance of a set of data points given a set of fitted coefficients. The deviance is given by:
	%
	%	D(y,mu) = 2\sum y_i ln (y_i / \mu_i) - y_i + \mu_i
	%
	%Usage:
	%	[y, tspks] = glmsim(model, data)
	%     
	%Input:
	%	model = a structure of fit coefficients from MLE_glmfit
	%	data = a structure of stimulus and spike history data from ./models
	%	maxspks = (optional, default = 0) If set to 1 then will enforce a maximum number of spikes per timestep equal
	%		to one spike per 2 milliseconds. 
	%   
	%Output:
	%	y = a vector of a simulated spike train given cursor data for each unit
	%  
	%Test code:
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	const = 'on';
	%	nK_sp = 50; 
	%	nK_pos = 10;
	%	dt_sp = 0.002;
	%	dt_pos = 0.05;
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%	model = MLE_glmfit(data, const);
	%	trains = glmsim(processed, model, data);

	if (nargin < 4) maxspks = 0; end
	bs = processed.binsize;
	%Force there to be a maximum number of spikes per time bin equal to binsize/0.001
	%meaning there is at most one spike per two milliseconds
	if maxspks == 1
		nMaxSpks = bs/0.002;
	end

	nU = size(data.X,1);
	N = size(data.X,2);
	nK_sp = length(data.sp_hist);
	y = zeros(N, nU);
	rho = zeros(N, nU);
	dev = zeros(nU);

	%indices of spike history filter
	sp_indices = data.sp_hist+1;
	%we're simulating new spikes, so we create a new spike history filter
	sp_hist = zeros(1,nK_sp);
	%time of spikes
	elecs = cell(1,nU);
	spikemuas = struct('times', elecs);
	for idx=1:nU
		spikemuas(idx).times = [0];    
	end
	%for each unit
	for i = 1:nU
		%extract all the unit's filters
		b_hat = model.b_hat(i,:);
		%and the spike history filter in particular
		k_sp = b_hat(sp_indices);
		%then set the spike history filter coefficients to zero, since we're generating new spikes and not using data.X's spike history
		b_hat(sp_indices) = 0;
		%compute the component of mu = e^eta that comes from the remaining filters used
		mu = glmval(b_hat', squeeze(data.X(i,:,:)), 'log');
		%then for each data point
		for j = 1:N
			%compute mu that incorporates the current spike history
			mu_sp = mu(j)*exp(sp_hist*k_sp');
			%then sample y ~ Pn(exp(eta)) to decide if we spike or not
			if maxspks == 1
				yij = min(poissrnd(mu_sp), nMaxSpks);
			else
				yij = poissrnd(mu_sp);
			end
			y(j,i) = yij;
			rho(j,i) = mu_sp;
			%If there's a spike, add the time(s) to spikemuas
			if (yij > 0)
				mu_sp;
				yij;
				if yij > 1e6
					display(['too many spikes per time bin: bin #: ' num2str(j) ' n. spikes: ' num2str(yij)])
				end
				sptime = bs*(j+0.005*(1:yij)');
				spikemuas(i).times = [spikemuas(i).times; sptime];
				%Only add this if there is a spike, otherwise set it to zero
				dev(i) = dev(i) + 2*yij*log(yij);
			end
			dev(i) = dev(i)-2*yij*log(mu_sp)-2*yij+2*mu_sp;
			%then update the spike history filter
			sp_hist = [sp_hist(2:end), yij];
		end
	end
	tspks = spikemuas;
