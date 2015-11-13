function [traindata, testdata] = makefolds(data, idx, fold)
	%Split data into training and test sets for cross validation of models
	%	Will split evenly into 'fold' number of components, and extract the 'idx'
	%	component of these folds as the test dataset, the rest being training data
	%
	%Usage:
	%	[traindata, testdata] = makefolds(data, idx, fold)
	%     
	%Input:
	%	data = data structure to split into test and train
	%	idx = current fold number to extract. Must be less than number of total folds
	%	fold = number of folds   
	%
	%Output:
	%	traindata = training data subset of data structure
	%	testdata = test data subset of data structure
	%  
	%Test code:
	%	pre = load('./testdata/test_preprocess_spline.mat');
	%	nK_sp = 50; 
	%	nK_pos = 1;
	%	dt_sp = 0.002;
	%	dt_pos = 0.2;
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%	idx = 1; fold = 10;
	%	%Split data into 9/10ths training, 1/10 test
	%	[train, test] = makefolds(data, idx, fold);

	N = size(data.y,2);
	assert(idx <= fold, 'Must set test data idx to be less than or equal to fold number');
	blocksize = floor(N/fold);
	preidx = 1:((idx-1)*blocksize);
	testidx = ((idx-1)*blocksize+1):(idx*blocksize);
	postidx = ((idx)*blocksize+1):N;
	trainidx = [preidx, postidx];
	traindata.y = data.y(:,trainidx);
	traindata.k = data.k;
	traindata.sp_hist = data.sp_hist;
	traindata.X = data.X(:,trainidx,:);
	traindata.torque = data.torque(trainidx,:);
	traindata.dtorque = data.dtorque(trainidx,:);
	traindata.ddtorque = data.ddtorque(trainidx,:);
	testdata.y = data.y(:,testidx);
	testdata.k = data.k;
	testdata.sp_hist = data.sp_hist;
	testdata.X = data.X(:,testidx,:);
	testdata.torque = data.torque(testidx,:);
	testdata.dtorque = data.dtorque(testidx,:);
	testdata.ddtorque = data.ddtorque(testidx,:);

