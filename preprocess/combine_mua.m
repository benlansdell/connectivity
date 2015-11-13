function mua_processed = combine_mua(sua_processed)
	%Combine single unit activity from a preprocess script to multi-unit activity of each electrode
	%
	%Usage:
	%	mua_processed = combine_mua(sua_processed)
	%
	%Input:
	%		sua_processed = output from preprocess script
	%		
	%Output:
	%		mua_processed is the same structure but with activity combined within the same electrode
	%
	%Test code:
	%		nevfile = './testdata/20130117SpankyUtah001.nev';
	%		binsize = 0.002;
	%		offset = 0.0;
	%		threshold = 5;
	%		fn_out = './worksheets/diagnostics/plots/test_spline_pre.eps';
	%		sua_processed = preprocess_spline(nevfile, binsize, threshold, offset, fn_out);
	%		mua = combine_mua(sua_processed);

    mua_processed = sua_processed;
    sua_electrodes = cellfun(@(x) floor(str2num(x)), sua_processed.unitnames);
    mua_electrodes = unique(sua_electrodes);
    mua_processed.unitnames = strtrim(cellstr(num2str(mua_electrodes'))');
    nE = length(mua_electrodes);
    N = size(sua_processed.binnedspikes,1);
	elecs = cell(1,nE);
    spikemuas = struct('times', elecs);
    mua_processed.binnedspikes = zeros(N, nE);
    mua_processed.rates = zeros(N, nE);
    %Find all units on the same electrode:
    for idx = 1:length(mua_electrodes)
    	electrode = mua_electrodes(idx);
    	indices = (sua_electrodes == electrode);
    	mua_processed.binnedspikes(:,idx) = sum(sua_processed.binnedspikes(:,indices),2);
    	mua_processed.rates(:,idx) = sum(sua_processed.rates(:,indices),2);
    	spikemuas(idx).times = [];
    	indicestocombine = find(indices);
    	for j = indicestocombine
    		spikemuas(idx).times = [spikemuas(idx).times; sua_processed.tspks(j).times];
    	end
    	spikemuas(idx).times = sort(spikemuas(idx).times);
    end
	mua_processed.tspks = spikemuas;