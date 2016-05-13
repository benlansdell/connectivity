function datout = truncateTrials(datin, timelimit)	
	%For GPFA analysis
	nT = length(datin)
	for idx = 1:nT
		datout(idx).trialId = datin(idx).trialId;
		datout(ntrials).quadrant = datin(idx).quadrant;
		datout(ntrials).oct = datin(idx).oct;
		nB = min(timelimit, size(datin(idx).spikes,2));
		datout(ntrials).spikes = datin(idx).spikes(:,1:nB);
	end
end