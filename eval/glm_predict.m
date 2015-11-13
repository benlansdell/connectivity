function [corrs, simTrain] = glm_predict(processed, data, model, fn_out, nTrains)
	%Predict spike trains for GLM given a model and cursor data. Compare the spike trains to the actual spike train generated
	%
	%Usage:
	%	y = glm_predict(processed, model, data, fn_out)
	%     
	%Input:
	%	processed = a structure output from a preprocess function
	%	data = a structure of stimulus and spike history data from ./models
	%	model = a structure of fit coefficients from MLE_glmfit
	%	fn_out = base file name to write plots to
	%	nTrains = (optional, default = 20) the number of trains to simulate for the same stimulus
	%   
	%Output:
	%	corrs = a vector of correlation coefficients between simulated firing rate, and actual firing rate
	%		 (a smoothed version of spike train)
	%	simTrain = average firing rate of all simulated trains. unsmoothed
	%	simTrains = all simulated trains
	%  
	%Test code:
	%	%Load test preprocessed data
	%	fn_out = './worksheets/09_12_2014/plots/test';
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	const = 'on';
	%	nK_sp = 50; 
	%	nK_pos = 10;
	%	dt_sp = 0.002;
	%	dt_pos = 0.05;
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%	model = MLE_glmfit(data, const);
	%	corrs = glm_predict(pre.processed, data, model, fn_out);

	nTrains = 20;
	%Number of trains to simulate monte-carlo style and estimate distribution of probability
	nTrainsMC = 100;
	nU = size(data.X,1);
	N = size(data.X,2);
	nK_sp = length(data.sp_hist);
	corrs = zeros(nU, 1);
	actualTrain = data.y';
	%For keeping track of the mean
	simTrains = zeros(nTrains, N, nU);
	%MC sims
	maxspks = 1;
	devianceMC = zeros(nTrainsMC, nU);
	devActual = deviance(model, data);
	smthfittedrates = zeros(nTrains, N, nU);
	%For keeping track of the square of activity (to compute the variance)
	%simTrain2 = zeros(N, nU);
	binsize = processed.binsize;

	%Make a Gaussian filter to smooth estimates
	sigma = 0.25;
	%sigma = 0.001;
	sigma = sigma/binsize;
	sz = sigma*3*3;
	x = linspace(-sz/2, sz/2, sz);
	gaussFilter = exp(-x.^2/(2*sigma^2));
	gaussFilter = gaussFilter/sum(gaussFilter);

	%Plot ten seconds worth of data
	t_i = 30;
	%t_f = 40;
	t_f = 50; 
	ii = 1:N;
	tt = ii*binsize;
	ii = ii(tt > t_i & tt < t_f);
	tt = tt(ii);

	for i = 1:nTrains
		i
		simt = glmsim(processed, model, data);
		simTrains(i,:,:) = simt;
		%simTrain2 = simTrain2 + simt.^2/(nTrains-1);
	end
	simTrain2 = std(simTrains, 1)';
	simTrain = mean(simTrains, 1)';
	%simTrain2 = sqrt(simTrain2-nTrains^2/(nTrains-1)^2*(simTrain).^2);

	display('Running MC simulations');
	for i = 1:nTrainsMC
		i
		[simt, tspks, rho, dev] = glmsim(processed, model, data, maxspks);
		devianceMC(i,:) = dev;
	end

	%for each unit
	for i = 1:nU
		for j = 1:nTrains
			smthfittedrates(j,:,i) = conv(simTrains(j,:,i), gaussFilter, 'same')/processed.binsize;
		end
		%Smooth this average train, and the original train
		smthfittedratesM = conv(simTrain(:,i), gaussFilter, 'same')/processed.binsize;
		smthrates = conv(actualTrain(:,i), gaussFilter, 'same')/processed.binsize;
		%Smooth the standard deviation, too
		smthfittedSD = conv(simTrain2(:,i), gaussFilter, 'same')/processed.binsize;
		%Compute correlation coefficient between the two trains
		actualSimCorr = corrcoef(smthfittedratesM, smthrates);
		corrs(i) = actualSimCorr(1,2);
		%Plot the simulated trains, the actual train, and the cursor data used to generate it
		clf;
		%Plot firing rate data
		%Plot rasters instead
		subplot(4,1,1)
		hold on
		spks = find(actualTrain(ii,i))*binsize+t_i;
		for k = 1:length(spks)
			plot([spks(k), spks(k)], [0 0.6], 'r', 'LineWidth', 0.1);
		end
		for j = 1:nTrains
			spks = find(simTrains(j,ii,i))*binsize+t_i;
			for k = 1:length(spks)
				plot([spks(k), spks(k)], [j; (j+0.6)], 'k', 'LineWidth', 0.1)
			end
		end
		ylim([0, nTrains+1])
		ylabel('trial')
		xlabel('time (s)')
		title(['Unit: ' processed.unitnames{i} ' correlation: ' num2str(corrs(i))]);

		subplot(4,1,2);
		hold on
		ymin = max(0,min(smthfittedratesM(ii)-smthfittedSD(ii))*1.2);
		ymax = max(smthfittedratesM(ii)+smthfittedSD(ii))*1.2
		%area(tt, smthfittedratesM(ii)+smthfittedSD(ii), ymin, 'FaceColor', [0.8 0.8 0.8])
		%area(tt, smthfittedratesM(ii)-smthfittedSD(ii), ymin, 'FaceColor', [1 1 1])
		plot(tt, smthfittedrates(:,ii,i), 'Color', [0.5 0.5 0.5])
		plot(tt, smthrates(ii), tt, smthfittedratesM(ii), 'LineWidth', 2)
		%plot(tt, data.y(idx,ii)/processed.binsize, tt, rho_hat(ii)/processed.binsize);
		%ylim([ymin, ymax])
		xlim([t_i, t_f])
		legend('Actual', 'GLM')
		xlabel('time (s)')
		ylabel('estimated firing rate (Hz)')

		%Plot cursor data
		subplot(4,1,3);
		plot(tt, data.torque(ii,1), tt, data.torque(ii,2))
		xlim([t_i, t_f])
		legend('RU', 'FE')
		xlabel('time (s)')		

		%Plot the distribution of deviances, along with the deviance of the actual simulation
		subplot(4,1,4);
		hist(devianceMC(devianceMC~=0));
		hold on
		plot(devActual, 0, '.r');
		pval = sum(devianceMC>devActual)/nTrainsMC;
		title(['p-value: ' num2str(pval)])
		xlabel('deviance');

		%Save plot
		saveplot(gcf, [fn_out '_unit_' processed.unitnames{i} '_pred.eps'], 'eps', [12 6]);
		%save fig
		saveas(gcf, [fn_out '_unit_' processed.unitnames{i} '_pred.fig'])
	end