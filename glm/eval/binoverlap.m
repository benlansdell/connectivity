function [angle_scores, angle_stdscores, angle_propabove, overlap_scores, overlap_stdscores, overlap_propabove, proptuned, proptunedsem, overlap_gcthresh, stdoverlap_gcthresh, nabove,	propbothru, propbothfe, propneitherferu, meantuningsizeangle, stdtuningsizeangle] = bindata(all_d, nbins)
	if nargin < 2
		nbins = 100;
	end
	%Bin them into bins...
	binned_overlap = cell(nbins,1);
	binned_angle = cell(nbins,1);
	%Set scores to zero
	binned_overlap = cell(nbins,1);
	binned_angle = cell(nbins,1);
	for idx = 1:size(all_d,1)
		%Figure out which bin to put score in...
		angle = all_d(idx,4);
		if angle < 0 
			angle = angle + 2*pi;
		end
		if angle > pi
			angle = angle - 2*pi;
		end
		overlap = cos(all_d(idx,4));
		score = all_d(idx,1);
		bin = ceil(nbins*(overlap+1)/(2));
		binned_overlap{bin} = [binned_overlap{bin},score];

		bin = ceil(nbins*(angle+pi)/(2*pi));
		binned_angle{bin} = [binned_angle{bin},score];
	end

	nK_sp = 6;
	sigthresh = 22.35; %Because: 1-chi2cdf(22.35, nK_sp) = 0.001;
	%sigthresh = 12.59; %Because: 1-chi2cdf(12.59, nK_sp) = 0.05;
	
	%Average scores, and variance of scores
	angle_scores = [];
	angle_stdscores = [];
	angle_propabove = [];
	overlap_scores = [];
	overlap_stdscores = [];
	overlap_propabove = [];
	for idx = 1:nbins
		angle_scores(idx) = mean(binned_angle{idx});
		angle_stdscores(idx) = std(binned_angle{idx});
		angle_propabove(idx) = sum(binned_angle{idx}>sigthresh)/length(binned_angle{idx});
		overlap_scores(idx) = mean(binned_overlap{idx});
		overlap_stdscores(idx) = std(binned_overlap{idx});
		overlap_propabove(idx) = sum(binned_overlap{idx}>sigthresh)/length(binned_overlap{idx});
	end

	gcthresholds = 0.1:5:200;
	anglethresh = 10/180*pi;
	proptuned = zeros(length(gcthresholds),1);
	diffangle = mod(all_d(:,2)-all_d(:,3), 2*pi);
	diffangle(diffangle < 0) = diffangle(diffangle < 0) + 2*pi;
	diffangle(diffangle > pi) = diffangle(diffangle > pi) - 2*pi;
	diffangle = abs(diffangle);
	
	quads1 = quadrants(all_d(:,2));
	quads2 = quadrants(all_d(:,3));

	bothfe = ((quads1 == 1) | (quads1 == 3)) & ((quads2 == 1) | (quads2 == 3));
	bothru = ((quads1 == 2) | (quads1 == 4)) & ((quads2 == 2) | (quads2 == 4));

	nabove = [];
	nsimtuned = [];
	stdproptuned = [];
	proptunedsem =[];
	overlap_gcthresh = [];
	stdoverlap_gcthresh = [];
	propbothru = [];
	propbothfe = [];
	propneitherferu = [];
	for idx = 1:length(gcthresholds)
		%Find those above threshold
		gcthresh = gcthresholds(idx);
		gcabove = all_d(:,1)>gcthresh;
		nabove(idx) = sum(gcabove);
		a1 = all_d(gcabove,2);
		a2 = all_d(gcabove,3);
		propbothfe(idx) = sum(bothfe(gcabove))/nabove(idx);
		propbothru(idx) = sum(bothru(gcabove))/nabove(idx);
		propneitherferu(idx) = 1-propbothfe(idx)-propbothru(idx);
		overlap_gcthresh(idx) = mean(cos(a1-a2));
		stdoverlap_gcthresh(idx) = std(cos(a1-a2))/sqrt(nabove(idx));
		%Find the prop of those whose difference in tuning angle is less than anglethresh
		nsimtuned(idx) = sum(diffangle(gcabove)<anglethresh);
		p = nsimtuned(idx)/nabove(idx);
		proptuned(idx) = p;
		stdproptuned(idx) = sqrt(p*(1-p));
		proptunedsem(idx) = stdproptuned(idx)/sqrt(nabove(idx));

	end

	if size(all_d, 2) == 5
		%Bin them into bins...
		tuningsizebinned = cell(nbins,1);
		for idx = 1:size(all_d,1)
			%Figure out which bin to put score in...
			tsize = all_d(idx,5);
			tangle = all_d(idx,2);
			if tangle < 0 
				tangle = tangle + 2*pi;
			end
			if angle > pi
				tangle = tangle - 2*pi;
			end
			bin = ceil(nbins*(tangle)/(2*pi));
			tuningsizebinned{bin} = [tuningsizebinned{bin},tsize];
		end
		meantuningsizeangle = [];
		stdtuningsizeangle = [];
		for idx = 1:nbins
			meantuningsizeangle(idx) = mean(tuningsizebinned{idx});
			stdntuningsizeangle(idx) = std(tuningsizebinned{idx});
		end
	else
		meantuningsizeangle = [];
		stdtuningsizeangle = [];
	end
end



