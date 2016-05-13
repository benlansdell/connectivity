function corrs = plot_predictions(model, data, processed, fn_out)
	%Plot use fitted GLM to predict firing rate of each unit based on filtered data
	%     
	%Input:
	%	model = data structure output by function in ./fitting (containing fitted coefficients)
	%	data = data structure output by ./models containing data used for fit
	%	processed = data structure output by ./preprocess containing processed raw data
	%	fn_out = (optional) base filename to write plots to for each unit. If not provided, just return correlations
	%
	%Output:
	%	corrs = correlation coefficient of each unit between smoothed actual and smoothed estimated firing rate
	%
	%Test code:
	%	const = 'on';
	%	fn_out = './worksheets/glm_ls/20130117SpankyUtah001';
	%	nK_sp = 100; 
	%	nK_pos = 0;
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
	%	%Fit model
	%	model = MLE_glmfit(data, const);
	%	%Plot predictions
	%	plot_predictions(model, data, pre.processed, fn_out);

	if (nargin < 4)
		plotflag = 0;
	else
		plotflag = 1;
	end
	nU = size(data.y,1); %number of units
	maxlags = 90/processed.binsize; %lag for x-correlation
	corrs = zeros(1,nU);

	%Make a Gaussian filter to smooth estimates
	sigma = 0.25;
	%sigma = 0.001;
	sigma = sigma/processed.binsize;
	sz = sigma*3*3;
	x = linspace(-sz/2, sz/2, sz);
	gaussFilter = exp(-x.^2/(2*sigma^2));
	gaussFilter = gaussFilter/sum(gaussFilter);

	%Plot ten seconds worth of data
	t_i = 30;
	%t_f = 33;
	t_f = 40; 
	ii = 1:size(data.y,2);
	tt = ii*processed.binsize;

	%For each unit, predicting firing rate, smooth and compare to actual smoothed firing rate
	for idx=1:nU 
		b_hat = model.b_hat(idx,:);
		%Uses the fit coefficients and the original input data to generate the ouput rho
		rho_hat = glmval(b_hat', squeeze(data.X(idx,:,:)), 'log');
		%Smooth estimates for plotting
		smthfittedrates = conv(rho_hat', gaussFilter, 'same')/processed.binsize;
		smthrates = conv(data.y(idx,:), gaussFilter, 'same')/processed.binsize;
		%Compute likelihood
		l = log_likelihood(rho_hat', data.y(idx,:));
		%Compute correlation between pred and actual smoothed firing rate
		pred_act_corr = corrcoef(smthfittedrates, smthrates);
		pred_act_corr = pred_act_corr(1,2);
		corrs(idx) = pred_act_corr;

		%Plot predictions
		if (plotflag ~= 0)
			clf;
			subplot(1,5,[1 2 3 4]);
			hold on
			plot(tt, smthrates(ii), tt, smthfittedrates(ii))
			%plot(tt, data.y(idx,ii)/processed.binsize, tt, rho_hat(ii)/processed.binsize);
			xlim([t_i, t_f])
			legend('Actual', 'GLM')
			xlabel('time (s)')
			ylabel('estimated firing rate (Hz)')
			title(['Unit: ' processed.unitnames{idx} '. Log-likelihood: ' num2str(l) '. Correlation: ' num2str(pred_act_corr)]);
	
			subplot(1,5,5);
			cc = xcorr(smthrates,smthfittedrates, maxlags, 'coeff');
			tt2 = linspace(-maxlags, maxlags, length(cc))*processed.binsize;
			plot(tt2,cc);
			xlabel('time (s)')
			title('cross-correlation')
	
			%save eps
			saveplot(gcf, [fn_out '_unit_' processed.unitnames{idx} '_fit.eps'], 'eps', [15,3.5]);  
			%save fig
			saveas(gcf, [fn_out '_unit_' processed.unitnames{idx} '_fit.fig'])
		end
	end
end