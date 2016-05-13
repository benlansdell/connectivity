function [processed1, processed2] = split_recording(processed, split, duration)
	%Split all channels into two parts
	%
	%Usage:
	%	[processed1, processed2] = split_recording(preprocessed, split, duration)
	%
	%Input:
	%		processed = output from preprocess script
    %       split = time point at which to split recording in two
    %       duration = total duration in seconds
	%		
	%Output:
	%		processed1/2 = the same structure but split
	%
	%Test code:
	%		nevfile = './testdata/20130117SpankyUtah001.nev';
	%		binsize = 0.002;
	%		offset = 0.0;
	%		threshold = 5;
    %       split = 300;
    %       duration = 360;
	%		processed = preprocess_spline(nevfile, binsize, threshold, offset);
    %       [p1 p2] = split_recording(processed, split, duration);

    short_processed = processed;
    nU = length(short_processed.unitnames);
    %Find number of bins to limit to
    nB = duration/processed.binsize;
    %Truncate processed structures
    short_processed.binnedspikes = short_processed.binnedspikes(1:nB,:);
    short_processed.rates = short_processed.rates(1:nB,:);
    if isfield(short_processed, 'torque')
        short_processed.torque = short_processed.torque(1:nB,:);
        short_processed.dtorque = short_processed.dtorque(1:nB,:);
        short_processed.ddtorque = short_processed.ddtorque(1:nB,:);
    end
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

    %Split further into components
    nB = split/processed.binsize;
    processed1 = short_processed;
    processed2 = short_processed;
    %Truncate processed structures
    processed1.binnedspikes = processed1.binnedspikes(1:nB,:);
    processed1.rates = processed1.rates(1:nB,:);
    if isfield(processed1, 'torque')
        processed1.torque = processed1.torque(1:nB,:);
        processed1.dtorque = processed1.dtorque(1:nB,:);
        processed1.ddtorque = processed1.ddtorque(1:nB,:);
    end
    %If cursor exists
    if isfield(processed1, 'dcursor')
        processed1.dcursor = processed1.ddcursor(1:nB,:);
        processed1.ddcursor = processed1.ddcursor(1:nB,:);
    end
    %If cursor exists
    if isfield(processed1, 'cursor')
        processed1.cursor = processed1.cursor(1:nB,:);
    end
    %If target info exists
    if isfield(processed1, 'target')
        processed1.target = processed1.target(1:nB,:);
    end
    %Cursor from trial structure
    if isfield(processed1, 'cursor_trial')
        processed1.cursor_trial = processed1.cursor_trial(1:nB,:);
    end
    %Remove all spike times past time point
    for idx = 1:length(nU)
        ltdur = (processed1.tspks(idx).times < split);
        processed1.tspks(idx).times = processed1.tspks(idx).times(ltdur);
    end

    %Truncate processed structures
    processed2.binnedspikes = processed2.binnedspikes(nB:end,:);
    processed2.rates = processed2.rates(nB:end,:);
    if isfield(processed2, 'torque')
        processed2.torque = processed2.torque(nB:end,:);
        processed2.dtorque = processed2.dtorque(nB:end,:);
        processed2.ddtorque = processed2.ddtorque(nB:end,:);
    end
    %If cursor exists
    if isfield(processed2, 'dcursor')
        processed2.dcursor = processed2.ddcursor(nB:end,:);
        processed2.ddcursor = processed2.ddcursor(nB:end,:);
    end
    %If cursor exists
    if isfield(processed2, 'cursor')
        processed2.cursor = processed2.cursor(nB:end,:);
    end
    %If target info exists
    if isfield(processed2, 'target')
        processed2.target = processed2.target(nB:end,:);
    end
    %Cursor from trial structure
    if isfield(processed2, 'cursor_trial')
        processed2.cursor_trial = processed2.cursor_trial(nB:end,:);
    end
    %Remove all spike times past time point
    for idx = 1:length(nU)
        ltdur = (processed2.tspks(idx).times > split) & (processed2.tspks(idx).times < duration);
        processed2.tspks(idx).times = processed2.tspks(idx).times(ltdur);
    end
end