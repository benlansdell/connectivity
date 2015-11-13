function short_processed = truncate_recording(processed, duration)
	%Truncate all channels to be no longer than a specified duration (in seconds)
	%
	%Usage:
	%	processed = truncate_recording(preprocessed, duration)
	%
	%Input:
	%		processed = output from preprocess script
    %       duration = duration in seconds
	%		
	%Output:
	%		processed = the same structure but truncated
	%
	%Test code:
	%		nevfile = './testdata/20130117SpankyUtah001.nev';
	%		binsize = 0.002;
	%		offset = 0.0;
	%		threshold = 5;
    %       duration = 360;
	%		processed = preprocess_spline(nevfile, binsize, threshold, offset);
    %       short_processed = truncate_recording(processed, duration);

    short_processed = processed;
    nU = length(short_processed.unitnames);
    %Find number of bins to limit to
    nB = duration/processed.binsize;
    %Truncate processed structures
    short_processed.binnedspikes = short_processed.binnedspikes(1:nB,:);
    short_processed.rates = short_processed.rates(1:nB,:);
    short_processed.torque = short_processed.torque(1:nB,:);
    short_processed.dtorque = short_processed.dtorque(1:nB,:);
    short_processed.ddtorque = short_processed.ddtorque(1:nB,:);
    %If cursor exists
    if isfield(short_processed, 'dcursor')
        short_processed.dcursor = short_processed.ddcursor(1:nB,:);
        short_processed.ddcursor = short_processed.ddcursor(1:nB,:);
    end
    %If cursor exists
    if isfield(short_processed, 'cursor')
        short_processed.cursor = short_processed.cursor(1:nB,:);
    end
    %If target info exists
    if isfield(short_processed, 'target')
        short_processed.target = short_processed.target(1:nB,:);
    end
    %Cursor from trial structure
    if isfield(short_processed, 'cursor_trial')
        short_processed.cursor_trial = short_processed.cursor_trial(1:nB,:);
    end
    %Remove all spike times past time point
    for idx = 1:length(nU)
        ltdur = (short_processed.tspks(idx).times < duration);
        short_processed.tspks(idx).times = short_processed.tspks(idx).times(ltdur);
    end
end