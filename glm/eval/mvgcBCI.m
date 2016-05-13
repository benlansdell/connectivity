function results = mvgcBCI(X, tstat, nK, alpha, momax, acmaxlags)

	if (nargin < 2) tstat = ''; end
	if (nargin < 3) nK = []; end
	if (nargin < 4) alpha = 0.05; end
	if (nargin < 5) momax = 20; end
	if (nargin < 6) acmaxlags = 1000; end
	
	demean = true;
	%%%%%%%
	%Setup%
	%%%%%%%

	if length(size(X))==3
		ntrials   = size(X,3);     % number of trials
		nobs      = size(X,2);   % number of observations per trial
		nvars = size(X,1);
	else
		ntrials = 1;
		nobs      = size(X,2);   % number of observations per trial
		nvars = size(X,1);
	end
	regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
	icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

	morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
	%momax     = 20;     % maximum model order for model order estimation

	%acmaxlags = momax;
	%acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

	%tstat     = 'chi2';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
	%alpha     = 0.05;   % significance level for significance test
	mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

	fres      = [];     % frequency resolution (empty for automatic calculation)	
	seed      = 0;      % random seed (0 for unseeded)


	%%%%%%%%%%%%%%%%%%%%%%
	%Estimate model order%
	%%%%%%%%%%%%%%%%%%%%%%

	if length(nK) == 0
		ptic('\n*** tsdata_to_infocrit\n');
		[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
		ptoc('*** tsdata_to_infocrit took ');
		
		results.AIC = AIC;
		results.BIC = BIC;
		results.moAIC = moAIC;
		results.moBIC = moBIC;

		fprintf('\nbest model order (AIC) = %d\n',moAIC);
		fprintf('best model order (BIC) = %d\n',moBIC);
		
		% Select model order.
		if strcmpi(morder,'actual')
		    morder = amo;
		    fprintf('\nusing actual model order = %d\n',morder);
		elseif strcmpi(morder,'AIC')
		    morder = moAIC;
		    fprintf('\nusing AIC best model order = %d\n',morder);
		elseif strcmpi(morder,'BIC')
		    morder = moBIC;
		    fprintf('\nusing BIC best model order = %d\n',morder);
		else
		    fprintf('\nusing specified model order = %d\n',morder);
		end
	else
		morder = nK;
	end

	results.morder = morder;

	%%%%%%%%%%%%%%%%%%%%%%
	%Fit VAR coefficients%
	%%%%%%%%%%%%%%%%%%%%%%

	ptic('\n*** tsdata_to_var... ');
	if (demean == true)
		display('Demeaning data')
		[A,SIG] = tsdata_to_var(X,morder,regmode);
	else
		display('Not demeaning data')
		[A,SIG] = tsdata_to_var_mean(X,morder,regmode);	
	end
	ptoc;

	results.var_A = A;
	results.var_SIG = SIG;
	% Check for failed regression
	assert(~isbad(A),'VAR estimation failed');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Auto covariance calculation%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	ptic('*** var_to_autocov... ');
	[G,info] = var_to_autocov(A,SIG,acmaxlags);	
	ptoc;
	results.G = G;
	results.info = info;
	var_info(info,true); % report results (and bail out on error)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Granger causality calculation%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Calculate time-domain pairwise-conditional causalities - this just requires
	% the autocovariance sequence.
	
	ptic('*** autocov_to_pwcgc... ');
	F = autocov_to_pwcgc(G);
	%[F, D] = autocov_to_pwcgc_chi2_ts(X, morder, regmode);
	ptoc;
	results.pwcgc = F;
	%results.pwcgc_chi2 = D;
	
	% Check for failed GC calculation
	assert(~isbad(F,false),'GC calculation failed');
	% Significance test using theoretical null distribution, adjusting for multiple
	% hypotheses.
	
	pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
	sig  = significance(pval,alpha,mhtc);

	%Compute chi2 pvals and significance

	%chi2pval = 1-chi2cdf(D, nK);
	%chi2sig = multiple_sig(chi2pval, alpha);
	%results.pwcgc_chi2pval = chi2pval;
	%results.pwcgc_chi2sig = chi2sig;

	results.pwcgc_pval = pval;
	results.pwcgc_sig = sig;

end