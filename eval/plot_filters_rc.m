function plot_filters_rc(model, data, processed, fn_out)
	%Plot filters of a fitted GLM model, along with other fit statistics.
	%Assumes that spike history filters are represented in a raised cosine basis
	%     
	%Input:
	%	model = data structure output by function in ./fitting (containing fitted coefficients)
	%	data = data structure output by ./models containing data used for fit
	%	processed = data structure output by ./preprocess containing processed raw data
	%	fn_out = base filename to write plots to for each unit
	%
	%Test code:
	%	const = 'on';
	%	fn_out = './worksheets/02_05_2015/20130117SpankyUtah001';
	%	nK_sp = 100; 
	%	nK_pos = 5;
	%	dt_pos = 0.2;
	%	dt_sp = 0.002;
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_55_1.mat');
	%	data = filters_sprc_pos(pre.processed, nK_sp, nK_pos);
	%	%Fit model
	%	model = MLE_glmfit(data, const);
	%	%Plot filters
	%	fn_out = './worksheets/glm_ls/20130117SpankyUtah001_rc';
	%	plot_filters_rc(model, data, pre.processed, fn_out);
	%	%Compare with usual basis
	%	datasp = filters_sp_pos(pre.processed, nK_sp, nK_pos);
	%	%Fit model
	%	modelsp = MLE_glmfit(datasp, const);
	%	%Plot filters
	%	fn_out = './worksheets/02_05_2015/20130117SpankyUtah001_sp';
	%	plot_filters(modelsp, datasp, pre.processed, fn_out);

	nU = size(data.y,1); %number of units
	k = data.k; %filter names and indices
	nK = size(k,1); %number of filters
	nP = 3; %number of things to plot about each fitted filter
	%For each unit, plot filters fit
	h = figure;
	for idx=1:nU 
		clf;
		b_hat = model.b_hat(idx,:);
		stats = model.stats{idx};
		%Extract spike history data
		rc_filt = b_hat(1+data.sp_hist);
		%Transform back to original basis
		sp_filt = data.spbasis*rc_filt';
		dt_filt = data.k{1,3};
		tt = (0:length(sp_filt)-1)*dt_filt*1000;
		%Find the location of the max pts of each column
		[maxspbasis,j] = find(data.spbasis == repmat(max(data.spbasis,[],1), size(data.spbasis,1),1));
		filtmaxs = maxspbasis*dt_filt*1000;
		if isfield(stats, 'se')	
			se = stats.se(k{1,2}+1)';
		else
			se = zeros(size(k{1,2}));
		end
		spse = data.spbasis*se';
		%Plot
		ymin = min(sp_filt-spse)*1.2;
		ymax = max(sp_filt+spse)*1.2;
		subplot(1,2,1)
		hold on
		area(tt, sp_filt+spse, ymin, 'FaceColor', [0.8 0.8 0.8])
		area(tt, sp_filt-spse, ymin, 'FaceColor', [1 1 1])
		plot(tt, sp_filt, 'LineWidth', 1, 'Color', [0 0 0]);
		plot(filtmaxs, zeros(size(filtmaxs)), 'go');
		spbasis = data.spbasis;
		for j=1:size(data.spbasis,2)
			spbasis(:,j) = spbasis(:,j).*rc_filt(j);
		end
		plot(tt, spbasis)
		if ischar(processed.unitnames)
			name = processed.unitnames;
		else
			name = processed.unitnames{idx};
		end
		title(['Unit: ' name ' no. of spikes: ' num2str(sum(processed.binnedspikes(:,idx)))]);
		xlabel('time (ms)');
		subplot(1,2,2)
		plot(tt, exp(sp_filt), 'LineWidth', 1, 'Color', [0 0 0]);
		title(['Unit: ' name ' no. of spikes: ' num2str(sum(processed.binnedspikes(:,idx)))]);
		xlabel('time (ms)');
		ylabel('gain')

		%save eps
		saveplot(gcf, [fn_out '_unit_' name '_filters.eps'], 'eps', [12,6]);	
	end
end