function status = csvwrite_heading(fn, M, titles)
	%csvwrite_heading	Write csv file for a matrix and a list of column headings.
	%
	% Usage:
	%			status = csvwrite_heading(fn, M, titles)
	%
	% Input:
	%			fn = output filename
	%			M = mxn matrix to output
	%			titles = {nx1} cell array containing strings of columne heading names
	%
	% Examples:
	%			titles = {'col a', 'col b', 'col c'};
	%			M = rand(3);
	%			csvwrite_heading('test.dat', M, titles);

	f_out = fopen(fn, 'w');
	%Write headings
	for idx = 1:length(titles)
		fprintf(f_out, [titles{idx} ',']);
	end
	fprintf(f_out, '\n');

	%Write the rest of the matrix
	for i = 1:size(M, 1)
		for j = 1:size(M,2)
			fprintf(f_out, [num2str(M(i,j)) ',']);
		end
		fprintf(f_out, '\n');
	end

	%Return with close status
	status = fclose(f_out);

end