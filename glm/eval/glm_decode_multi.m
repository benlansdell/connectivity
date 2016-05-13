function output = glm_decode_multi(processed, data, model, F, Q, mu, R, P, d, fn_out)
	%Predict cursor position using recursive Bayesian filter for GLM given a model 
	%	compares to cursor data.
	%
	%Usage:
	%	output = glm_decode_multi(processed, data, model, F, Q, mu, fn_out)
	%     
	%Input:
	%	processed = a structure output from a preprocess function
	%	data = a structure of stimulus and spike history data from ./models
	%	model = a structure of fit coefficients from MLE_glmfit
	%	F = cursor dynamics from AR model
	%	Q = covariance matrix from AR model
	%	mu = average cursor data from AR model
	%	R = order of AR model
	%	P = number of future points to predict AR model at (encoding model requires _future_ time points)
	%	d = resolution (in multiples of dt_sp) AR model is defined at
	%	fn_out = base file name to write plots to
	%   
	%Output:
	%	output = Predicted decoded cursor from spike trains
	%  
	%Test code:
	%	%Load test preprocessed data
	%	fn_out = './worksheets/11_17_2014/testdecode.eps';
	%	pre = load('./testdata/test_preprocess_spline_short24.mat');
	%	%pre = load('./testdata/test_preprocess_spline.mat');
	%	const = 'on';
	%	nK_sp = 50; 
	%	nK_pos = 4;
	%	dt_sp = 0.002;
	%	dt_pos = 0.2;
	%	R = 5;
	%	P = nK_pos;
	%	d = dt_pos/dt_sp;
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%	model = MLE_glmfit(data, const);
	%	[F, Q, mu] = fit_AR_LS_multi(data.torque, R, P, d);
	%	%load('./testdata/testglm.mat')
	%	decoded_torque = glm_decode_multi(pre.processed, data, model, F, Q, mu, R, P, d, fn_out);

	binsize = processed.binsize;
	output = zeros(size(data.torque));
	nU = size(data.X,1);
	N = size(data.X,2);
	epsilon = 1e-7;
	N_sub = 2000;
	%Only decode a short sample...
	N = 20000;

	%Make a Gaussian filter to smooth estimates
	sigma = 0.25;
	%sigma = 0.001;
	sigma = sigma/binsize;
	sz = sigma*3*3;
	x = linspace(-sz/2, sz/2, sz);
	gaussFilter = exp(-x.^2/(2*sigma^2));
	gaussFilter = gaussFilter/sum(gaussFilter);

	%Raise error if filter lengths are greater than 1... not implemented yet :(
	%if length(data.k{2,2}) > 1 | length(data.k{3,2}) > 1
	%	error('NotImplementedError: cursor filters must currently be of length 1')
	%end

	%GLM params
	%Cursor filters
	kx = model.b_hat(:,data.k{2,2}+1);
	ky = model.b_hat(:,data.k{3,2}+1);

	Q = sparse(Q + epsilon*eye(size(Q)));
	%Q(2*P+1:end,2*P+1:end) = eye(2*d*(R-1))*epsilon;

	%Our initial estimate
	%xk = mu';
	%The actual point
	xk = zeros(2*P+2*d*(R-1),1);
	xk(1:2:(2*P-1)) = data.torque(d*(R-1)+1+(0:d:(d*(P-1))),1);
	xk(2:2:(2*P)) = data.torque(d*(R-1)+1+(0:d:(d*(P-1))),2);
	xk(2*P+(1:2:(2*d*(R-1)-1))) = data.torque((d*(R-1)):-1:1,1);
	xk(2*P+(2:2:(2*d*(R-1)))) = data.torque((d*(R-1)):-1:1,2);
	Wk = Q;
	F = sparse(F);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%X = reshape([X_RU;X_FE],1,[]);
	%gradloglambda = [kx, ky];
	gradloglambda = sparse(zeros(nU,2*P+2*d*(R-1)));
	gradloglambda(:,1:(2*P)) = reshape([kx; ky], nU, []);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%timing and percentage indicator
	time_pu = 0;
	time_inv = 0;
	p = 0;
	%At each time step
	for idx = (d*(R-1)+1):N
	%for idx = (d*(R-1)+2):(d*(R-1)+20)
		if floor(idx*100/N) ~= p
			display([num2str(floor(idx*100/N)) '% decoding completed'])
			p = p + 1;
		end
		tic	
		%Do prediction step 
		xk_p = mu+F*xk;
		Wk_p = F*Wk*F'+Q;
		%Wk_p = f*Wk*f+Q;
		%Do update of covariance
		updateW = sparse(size(Wk,1),size(Wk,2));
		updateX = 0;
		for j = 1:nU
			%j=9;
			%Evaluate GLM for each unit, at given time step
			%replace in data.X the actual cursor data for what's been predicted
			x = squeeze(data.X(j,idx,:));
			RUidx = data.k{2,2};
			FEidx = data.k{3,2};
			x(RUidx) = xk_p(1:2:2*P-1);
			x(FEidx) = xk_p(2:2:2*P);
			b_hat = model.b_hat(j,:);
			lambda = glmval(b_hat', x', 'log');

			%Sum update contributions for each unit
			%updateW = updateW + lambda*gradloglambda(j,:)'*gradloglambda(j,:)*binsize;
			%updateX = updateX + gradloglambda(j,:)'*(data.y(j,idx)-lambda*binsize);
			updateW = updateW + lambda*gradloglambda(j,:)'*gradloglambda(j,:);
			updateX = updateX + gradloglambda(j,:)'*(data.y(j,idx)-lambda);
		end
		time_pu = time_pu + toc;
		%Don't let updateW get too big... or Wkinv will get too small...
		%if updateW > 1e10
		%	display('Very large W')
		%end
		%updateW = min(updateW, 1e10);

		tic
		Wkinv = inv(Wk_p)+updateW;
		Wk = inv(Wkinv);
		%Update position
		xk = xk_p+Wk*updateX;
		time_inv = time_inv + toc;
		%display(['At bin=' num2str(idx) ' cond(Wk)=' num2str(condest(Wk))]);
		display(['At bin=' num2str(idx)]);
		%Save result
		output(idx,:) = xk(1:2);
		%display(['idx=' num2str(idx) ', time=' num2str(idx*binsize) ', lambda=' num2str(lambda)])
		%xk
		%Wk 
		%pause
	end

	%Plot ten seconds worth of data
	t_i = 0;
	%t_f = 40;
	t_f = 40; 
	ii = 1:N;
	tt = ii*binsize;
	ii = ii(tt > t_i & tt < t_f);
	tt = tt(ii);

	%Plot cursor data
	subplot(2,1,1);
	plot(tt, output(ii,1), 'r--', tt, data.torque(ii,1), 'k')
	xlim([t_i, t_f])
	ylabel('Pred. RU')
	%Take a random sample... so it doesn't take forever...
	ii_sub = datasample(ii, N_sub, 'Replace', false);
	corrRU = corr(output(ii_sub,1), data.torque(ii_sub,1))
	corrFE = corr(output(ii_sub,2), data.torque(ii_sub,2))
	title(['Correlation RU:' num2str(corrRU) ' Correlation FE: ' num2str(corrFE)])
	subplot(2,1,2);
	plot(tt, output(ii,2), 'r--', tt, data.torque(ii,2), 'k')
	xlim([t_i, t_f])
	ylabel('Pred. FE')
	xlabel('time (s)')		
	
	saveplot(gcf, fn_out, 'eps', [9 6]);
	saveas(gcf, [fn_out])
end