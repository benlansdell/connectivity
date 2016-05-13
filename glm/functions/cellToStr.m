function colstr = cellToString(cellarray, quotationmark)
	if nargin < 2
		qm = '`';
	else 
		qm = quotationmark;
	end

	colstr = '';
	for idx = 1:length(cellarray)
		if idx > 1
			colstr = [colstr ',' qm cellarray{idx} qm];
		else
			colstr = [colstr qm cellarray{idx} qm];
		end
	end
end

