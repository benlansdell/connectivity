function plot_filters(model, data, processed, fn_out)
	%Plot filters of a fitted GLM model, along with other fit statistics
	%     
	%Input:
	%	model = data structure output by function in ./fitting (containing fitted coefficients)
	%	data = data structure output by ./models containing data used for fit
	%	processed = data structure output by ./preprocess containing processed raw data
	%	fn_out = base filename to write plots to for each unit
	%
	%Test code:
	%	const = 'on';
	%	fn_out = './worksheets/glm_ls/20130117SpankyUtah001';
	%	nK_sp = 100; 
	%	nK_pos = 100;
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
	%	%Fit model
	%	model = MLE_glmfit(data, const);
	%	%Plot filters
	%	plot_filters(model, data, pre.processed, fn_out);

	nU = size(data.y,1); %number of units
	k = data.k; %filter names and indices
	nK = size(k,1); %number of filters
	nP = 3; %number of things to plot about each fitted filter
	%For each unit, plot filters fit
	h = figure;
	for idx=1:nU 
		clf;
		b_hat = model.b_hat(idx,:);
		dev = model.dev{idx};
		stats = model.stats{idx};
		const = b_hat(1);
		for j = 1:nK
			%Extract data
			name = k{j,1};
			filt = b_hat(k{j,2}+1);
			%If available, find statistics of fit, otherwise set these to zero
			if isfield(stats, 'se')	
				se = stats.se(k{j,2}+1)';
				tstat = stats.t(k{j,2}+1);
				pval = stats.p(k{j,2}+1);
			else
				se = zeros(size(k{j,2}));
				tstat = zeros(size(k{j,2}));
				pval = zeros(size(k{j,2}));
			end
			dt_filt = k{j,3};
			%If filter length is zero skip this one
			if length(k{j,2}) < 1
				continue
			end

			%Plot filter plus/minus SE
			subplot(nP, nK+1, j)
			tt = (0:length(filt)-1)*dt_filt*1000;
			hold on
			ymin = min(filt-se)*1.2;
			ymax = max(filt+se)*1.2;
			area(tt, filt+se, ymin, 'FaceColor', [0.8 0.8 0.8])
			area(tt, filt-se, ymin, 'FaceColor', [1 1 1])
			plot(tt, filt);
			ylim([ymin ymax]);
			if length(tt) > 1
				xlim([min(tt) max(tt)]);
			end
			title(name);
			%xlabel('time (ms)');
			%Plot tstats
			subplot(nP, nK+1, (nK+1)+j)
			plot(tt, tstat);
			if (j == 1)
				ylabel('t statistic');
			end
			%xlabel('time (ms)');
			%p-values
			subplot(nP, nK+1, (nK+1)*2+j)
			hold on
			above = pval > 0.05;
			below = pval <= 0.05;
			plot(tt(above), log(pval(above)), '.b');
			plot(tt(below), log(pval(below)), '.r');
			plot(tt, log(0.05)*ones(length(tt),1), 'k')
			if (j == 1)
				ylabel('log(p-val)');
			end
			xlabel('time (ms)');
		end

		if ischar(processed.unitnames)
			name = processed.unitnames;
		else
			name = processed.unitnames{idx};
		end

		%Plot information about each subplot
		subplot(nP, nK+1, (nK+1))
		if isfield(stats, 'dfe')
			str1(1) = {['Unit: ' name]};
			str1(2) = {['Deviance: ' num2str(dev)]};
			str1(3) = {['Degrees of freedom: ' num2str(stats.dfe)]};
			str1(4) = {['Estimated dispersion: ' num2str(stats.sfit)]};
			str1(5) = {['Binsize: ' num2str(processed.binsize)]};
			str1(6) = {['Seconds of training: ' num2str(size(data.y,2)*processed.binsize)]};
			str1(7) = {['Number of spikes: ' num2str(sum(data.y(idx,:)))]};
			text(0.1,0.8,str1)
			axis off
		end
		subplot(nP, nK+1, (nK+1)*2)
		str2(1) = {'T-test'};
		text(0.1,0.8,str2)
		axis off
		subplot(nP, nK+1, (nK+1)*3)
		str3(1) = {'Test for significance of each predictor'};
		text(0.1,0.8,str3)
		axis off

		%save eps
		saveplot(gcf, [fn_out '_unit_' name '_filters.eps'], 'eps', [12,6]);	
		%save fig
		saveas(gcf, [fn_out '_unit_' name '_filters.fig'])
	end