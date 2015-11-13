function sig = multiple_sig(pvals, alpha, method)
	%Performs multiple hypothesis test using either
	%		- Bonferroni correction at alpha
	%		- Benjamini-Hochberg FDR control at alpha
	%		- basic thresholding (no multiple hypothesis correction)
	%
	%Usage:
	%		sig = multiple_sig(pvals, method, alpha)
	%
	%Input:
	%		pvals = vector of pvalues to be multiple hypothesis test corrected
	%		alpha = level at which to control false discoveries
	%		method = (optional, default = 'fdr') one of 'fdr', 'bonferroni', 'none'
	%
	%Test code:
	%		pvals = rand(100,1);
	%		pvals = [pvals; 1e-9];
	%		sig = multiple_sig(pvals, 0.05);

	if (nargin < 3) method = 'fdr'; end

	sig = nan;
	if strcmp(method, 'bonferroni')
		m = numel(pvals);
		thresh = alpha/m;
		sig = pvals < thresh;
	elseif strcmp(method, 'none')
		sig = pvals < alpha;
	else
		if ~strcmp(method, 'fdr')
			display(['multiple_sig: method ' char(method) ' not recognized. Specify one of fdr or bonferroni or none.'])
			display(['Continuing with fdr method.'])
		end
		[h crit_p adj_p]=fdr_bh(pvals,alpha);
		sig = h;		
	end
end
