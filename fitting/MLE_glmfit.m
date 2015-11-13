function model = MLE_glmfit(data, const)
	%Fit GLM to spike data from blackrock recording file for each unit above a specified threshold
	%     
	%Input:
	%	data = covariate data output structure from any function in ./models
	%	const = (optional, default = 'on') whether to fit a constant term to the model or not, 'on' or 'off'
	%   
	%Output:
	%	model is a structure containing the following fields:
	%		b_hat = [nU x (nK + 1)] array with spikes from all channels binned according to binsize. nB = no. bins, nU = no. units.
	%			Note: if a constant term is not fit, a column of zeros is appended to b_hat to make dimensions consistent
	%		dev = [nU x 1] cell array listing deviance of each unit's fit
	%		stats = [nU x 1] cell array listing fitting statistics output from glmfit
	%		converged = [nU x 1] array listing 1 if the IRLS converged within iteration limit, 0 if not
	%		conditioned = [nU x 1] array listing 1 if the IRLS did not issue an ill-conditioned warning, 0 if it did
	%
	%Test code:
	%	const = 'on';
	%	nK_sp = 100; 
	%	nK_pos = 100;
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
	%	model = MLE_glmfit(data, const);

	if (nargin < 2) const = 'on'; end

	nU = size(data.y,1);
	nK = size(data.X,3);
	if strcmp(const, 'on')
		model.b_hat = zeros(nU, nK+1);
		model.mask = zeros(nU, nK+1);
	else
		model.b_hat = zeros(nU, nK);
		model.mask = zeros(nU, nK);
	end
	model.dev = cell(nU,1);
	model.stats = cell(nU,1);
	model.converged = ones(nU,1);
	model.conditioned = ones(nU,1);
	cutoff = 1e-10;
	%For each unit, fit a GLM to the torque data
	display(['Fitting GLM by MLE with IRLS. Fitting ' num2str(nU) ' units.'])
	for idx=1:nU 
		display(['Fitting unit ' num2str(idx)])
		%Mask columns that don't vary... they cannot be estimated.
		d = squeeze(data.X(idx,:,:));
		mask = (std(d) > cutoff);
		if strcmp(const, 'off')
			m = mask;
		else
			m = [1==1 mask];
		end
		model.mask(idx,:) = m;
		[b, dev, stats] = glmfit(d(:,mask),data.y(idx,:),'poisson', 'constant', const);
		%Catch if a warning was raised about badly conditioned matrix
		[warn, warnid] = lastwarn;
		if ~strcmp(warn, '')
	   		switch warnid
        	case 'stats:glmfit:IterationLimit'
        		model.converged(idx) = 0;
        	case 'stats:glmfit:BadScaling'
        		model.conditioned(idx) = 0;
       		end
	    end
	    lastwarn('')
		%Extract filters fitted...
		model.b_hat(idx,m) = b;	
		se = zeros(size(model.b_hat,2));
		se(idx,m) = stats.se;
		stats.se = se;
		model.dev{idx} = dev;
		%Remove residual components since these take up a lot of memory
		model.N = size(stats.resid,1);
		stats = rmfield(stats, {'resid', 'residp', 'residd', 'resida', 'wts'});
		model.stats{idx} = stats;
	end
	model.logli = ll(model, data, 'poisson');
	if ~strcmp(const, 'on')
		model.b_hat = [zeros(nU, 1), model.b_hat]
	end
	display('Done')
