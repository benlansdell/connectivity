function [GCdev, GCpval, GCsig, fulldeviances] = granger(processed, data, fn_out, pval)
	%Computes Granger-causality matrices for set of spike and cursor trajectory data
	%
	%Usage:
	%	[GCdev, GCpval, GCsig] = granger(processed, data, fn_out, pval);
	%     
	%Input:
	%	processed = a structure
	%	data = a structure of stimulus and spike history data from ./models directory
	%	fn_out = filename to output plots
	%	pval = (optional, default=0.001) p-value at which to perform statistical significance test 
	%   
	%Output:
	%	GCdev = [nU x nU] matrix where i,j element is the change in deviance when unit i is excluded
	%		from a model of unit j. (i.e. unit i's effect on unit j)
	%	GVpval = [nU x nU] matrix listing p-value of each GCdev value
	%	GVsig = [nU x nU] matrix listing significance of GCdev
	%  
	%Test code:
	%	%Load test preprocessed data
	%	const = 'on';
	%	pval = 0.001;
	%	nK_sp = 6; 
	%	nK_pos = 6;
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_60hz_short24.mat');
	%	fn_out = './worksheets/12_2_2014/plots/testgranger.eps';
	%	processed = pre.processed;
	%	data = filters_sp_pos_network(processed, nK_sp, nK_pos);
	%	[GCdev, GCpval, GCsig] = granger(processed, data, fn_out, pval);

	if (nargin < 4) pval = 0.001; end
	nU = size(data.y,1);
	nK = size(data.X,2);
	nKL = size(data.k,1);
	nK_sp = size(data.k{1,2},1);

	const = 'on';
	%Estimate unperturbed model
	model = MLE_glmfit_network(data, const);
	%Add deviance to GCdev
	devs = zeros(1,nU);
	for i = 1:nU
		devs(i) = -model.dev{i};
	end
	fulldeviances = -devs;
	GCdev = repmat(devs,nU,1);

	%Modify data matrix in turn for each unit
	for i = 1:nU
		%Modify data matrix so that activity for unit i is left out
		notidx = data.k{i,2};
		startidx = notidx(1)-1;
		endidx = notidx(end)+1;
		idx = [1:startidx, endidx:nK];
		data_noti = data;
		data_noti.X = data_noti.X(:,idx);
		%Should really also update the data.k information to match what's going on, 
		%but this isn't used by the MLE_glmfit function, and nothing else downstream
		%is done with the model, currently, so this is fine.
		%Estimate models
		model = MLE_glmfit_network(data_noti, const);
		devs = zeros(1,nU);
		for j = 1:nU
			devs(j) = model.dev{j};
		end
		%Record change in deviance
		GCdev(i,:) = GCdev(i,:)+devs;
	end

	%Compute p-values of each GC measure
	GCpval = 1-chi2cdf(GCdev, nK_sp);
	%Threshold to compute result of significance test
	%pval = 0.05;
	%GCsig = GCpval < pval;
	GCsig = multiple_sig(GCpval, pval);

	for i = 1:nU
		GCdev(i,i) = 0;
		GCsig(i,i) = 0;
	end

	if isstr(fn_out)
		%Plot results
		clf
		subplot(3,1,1)
		colormap(bone);
		imagesc(GCdev)
		title('Change in deviance')
		ylabel('Unit')
		xlabel('Unit')
		set(gca,'XTick',1:nU);
		set(gca,'YTick',1:nU);
		set(gca,'XTickLabel',processed.unitnames);
		set(gca,'YTickLabel',processed.unitnames);
		rotateXLabels(gca, 90);
		colorbar
		subplot(3,1,2)
		imagesc(GCpval)
		title(['p-value'])
		ylabel('Unit')
		xlabel('Unit')
		set(gca,'XTick',1:nU);
		set(gca,'YTick',1:nU);
		set(gca,'XTickLabel',processed.unitnames);
		set(gca,'YTickLabel',processed.unitnames);
		%rotate x labels
		rotateXLabels(gca, 90);
		colorbar
		subplot(3,1,3)
		imagesc(GCsig)
		title(['Significance test @ p<' num2str(pval)])
		ylabel('Unit')
		xlabel('Unit')
		set(gca,'XTick',1:nU);
		set(gca,'YTick',1:nU);
		set(gca,'XTickLabel',processed.unitnames);
		set(gca,'YTickLabel',processed.unitnames);
		rotateXLabels(gca, 90);
		colorbar
		saveplot(gcf, fn_out, 'eps', [6 18])
	end
	%Cluser the results and plot the clustered matrices
	%fn_out_c = [fn_out '.cluster'];
	%granger_cluster(GCdev, GCpval, GCsig, processed.unitnames, fn_out_c);