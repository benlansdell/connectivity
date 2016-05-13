function data = filters_sprc_pos_network(processed, nK_sp, nK_pos, dt_sp, dt_pos)
	%Prepare spike and torque data for GLM which includes spike history and cursor position (x_1, x_2) filters:
	%Spike history filter is saved as a rasied cosine function coefficients
	%
	%	y(i) ~ Pn(g(eta_i))
	%
	%where
	%
	%	eta_i = \sum y(i-j) k_sp(i) + \sum x_1(i+j) k_1(j) + \sum x_2(i+j) k_2(j)
	%
	%Usage:
	%	data = filters_sprc_pos_network(processed, unit, nK_sp, nK_pos, dt_sp, dt_pos)
	%     
	%Input:
	%	processed = structure output from one of the preprocess functions.
	%	unit = index of unit to compute data matrix for. Large matrix sizes suggest we compute one unit at 
	%		a time
	%	nK_sp = number of timebins used for spike history filter for all units
	%	nK_pos = number of timebins used for cursor trajectory filters (on in x and y axis)
	%	dt_sp = (optional, default = binsize in processed structure) step size of spike history filters
	%		in seconds. Must be a multiple of the data's binsize.
	%	dt_pos = (optional, default = binsize in processed structure) step size of position filter in
	%		seconds. Must be a multiple of the data's binsize
	%   
	%Output:
	%	data is a structure containing the following fields:
	%		y = [nU x nB] array where y_ij is the number of spikes at time bin j for unit i.
	%		X = [nB x nK] array where X_ijk is the value of covariate k, at time bin j, for unit i
	%			Note: nK = nU*nK_sp + 2*nK_pos
	%		k = Names of each filter, a [n x 2] cell array in which each row is of the form ['filter j name', [idxj1 idxj2 ...]]
	%			Note: The second column lists indices in 1:nK to which the label applies
	%		torque = torque data trimmed in the same way X and y are. 
	%			Note: truncated at start and end because spike and cursor trajectory are not defined for first 
	%			and last nK_sp and nK_pos timebins respectively.
	%		dtorque = trimmed dtorque
	%		ddtorque = trimeed ddtorque
	%  
	%Test code:
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_60hz_short24.mat');
	%	nK_sp = 100; 
	%	nK_pos = 5;
	%	dt_sp = 0.002;
	%	dt_pos = 0.2;
	%	data = filters_sprc_pos_network(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);

	if (nargin < 5) dt_sp = processed.binsize; end
	if (nargin < 6) dt_pos = processed.binsize; end

	%Check dt's specified are valid
	assert(rem(dt_sp,processed.binsize)==0, 'Invalid dt_sp. Must be a multiple of binsize');
	assert(rem(dt_pos,processed.binsize)==0, 'Invalid dt_pos. Must be a multiple of binsize');
	steps_sp = dt_sp/processed.binsize;
	steps_pos = dt_pos/processed.binsize;

	nU = size(processed.binnedspikes,2);
	nB = size(processed.binnedspikes,1);
	nK = nU*nK_sp + 2*nK_pos;

	T = nK_sp*dt_sp;
	[rcbasis, spbasis, nK_rc] = makeRCBasis(dt_sp, T);
	nK = nU*nK_rc + 2*nK_pos;

	data.X = zeros(nB, nK);
	data.k = cell(nU+2,3);
	for i = 1:nU
		data.k{i,1} = ['spike unit' num2str(i)]; 
		data.k{i,2} = (i-1)*nK_rc + (1:nK_rc);
		data.k{i,3} = dt_sp;
	end
	data.k{nU+1,1} = 'RU pos'; 
	data.k{nU+1,2} = (1:nK_pos) + nU*nK_rc;
	data.k{nU+1,3} = dt_pos;
	data.k{nU+2,1} = 'FE pos'; 
	data.k{nU+2,2} = (1:nK_pos) + nU*nK_rc + nK_pos;
	data.k{nU+2,3} = dt_pos;
	%Record specifically which indices are spike history indices for model simulation
	%data.sp_hist = data.k{1,2};

	%For each unit, add data to X array
	%Make stimulus vector at each timebin
	for j = (nK_sp*steps_sp+1):(nB-nK_pos*steps_pos)
		shist = zeros(nK_rc*nU,1);
		for idx=1:nU 
			%(past) spike history
			shist(((idx-1)*nK_rc+1):(idx*nK_rc)) = project_rc(processed.binnedspikes(j-nK_sp*steps_sp:steps_sp:j-steps_sp, idx), rcbasis);
		end
		%(future) torque trajectory
		torqueRU = processed.torque(j:steps_pos:(j+(nK_pos-1)*steps_pos),1);
		torqueFE = processed.torque(j:steps_pos:(j+(nK_pos-1)*steps_pos),2);
		%Form stim vector
		data.X(j,:) = [shist' torqueRU' torqueFE'];
	end
	%Truncate to exclude start and end of recording where spike history 
	%and cursor trajectory aren't well defined
	data.X = data.X((nK_sp*steps_sp+1):(nB-nK_pos*steps_pos),:); %(nkt+1:end-nkt,:);
	data.y = processed.binnedspikes((nK_sp*steps_sp+1):(nB-nK_pos*steps_pos), :)';
	%Truncate other data for comparison, too
	data.torque = processed.torque((nK_sp*steps_sp+1):(nB-nK_pos*steps_pos),:); 
	data.dtorque = processed.dtorque((nK_sp*steps_sp+1):(nB-nK_pos*steps_pos),:);
	data.ddtorque = processed.ddtorque((nK_sp*steps_sp+1):(nB-nK_pos*steps_pos),:);

	data.rcbasis = rcbasis;
	data.spbasis = spbasis;
end

function sphistory_rc = project_rc(sphistory, rcbasis)
	sphistory_rc = rcbasis*sphistory;
end

function [rcbasis, spbasis, nK_rc] = makeRCBasis(dt, T)
	%Define basis of raised cosine functions
	%Hard-wired log-scale
	a = 15;
	nTotal = floor(T/dt);
	%Create basis function function
	B = @(t, a, psi, phi) iif(a*log(t-psi)>phi-pi & a*log(t-psi)<phi+pi, 0.5*(cos((a*log(t-psi)-phi))+1), ...
		true, zeros(size(t)));
	tt = 0:dt:(T-dt);
	nT = length(tt);
	%Compute phi0 and psi for desired basis
	phi1 = a*log(dt*(1-exp(-pi/a))^(-1));
	phi0 = phi1-pi;
	psi = -exp(phi0/a);
	%Compute each function for phi_i
	phi = phi0;
	%Compute matrix of basis vectors
	rc = zeros(nT, nTotal);
	for j = 1:nTotal
		for i = 1:nT
			rc(i, j) = B(tt(i), a, psi, phi);
		end
		%Normalize each column
		if norm(rc(:,j))>0
			rc(:,j) = rc(:,j)/norm(rc(:,j));
		end
		phi = phi + pi;
	end
	%Truncate to just the non-zero columns
	nK_rc = rank(rc);
	rc = rc(:,1:nK_rc);
	%Flip so time near spike is best resolved
	rc = flipud(rc);
	%Compute pseudo inverse of this matrix
	rcbasis = inv(rc'*rc)*rc';
	spbasis = rc;
end