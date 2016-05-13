function saveplot(h, filename, format, dimensions, journal)
%saveplot	Saves the specified plot.
%
% Usage:
%			saveplot(h, filename, format, dimensions, fontsize)
%
% Input:
%			h = figure object. Can just use 'gcf'
%			filename = output filename
%			format = one of 'eps', 'png', 'jpg', 'pdf' or 'svg' (optional, default 'eps')
%			dimensions = (optional, default = [6, 4]) width and height of figure in inches OR
%					a scale factor to multiply the journal specific dims described below
%			journal = (optional) if specified as one of 'plos1', 'plos2', 'plos1.5' or 'els1', 'els2', 'els1.5'
%					will set width to that journals specified one column, two column, 1.5 column width (inches)
%					which can then be multiplied by the scale factor in dimensions if this plot is part
%					of a subplot -- otherwise just set dimensions to 1.
%
% Examples:
%			plot(0:10, (0:10).^2);
%			xlabel('time');
%			saveplot(gcf, './test.eps', 'eps', 1, 'plos1');

%Note that matlab saveas and open can be used to store figures for editing later...

	if (nargin < 3)
		format = 'eps';
		dimensions = [6 4];
	elseif nargin < 4
		dimensions = [6 4];
	elseif nargin == 4
		if length(dimensions) == 1
			throw(MException('Argument:Error', 'If specifying a scale dimension, must also specify journal. See help saveplot'));
		end
	elseif nargin == 5
		%If we're specifying the journal...get the axes handle
		ax = findall(h, 'type', 'axes');
		if length(dimensions) ~= 1
			throw(MException('Argument:Error', 'If specifying a journal, must specify scale factor. See help saveplot'));
		end
		%Measurements in inches
		if strcmp(journal, 'plos1')
			width = 3.27;
		elseif strcmp(journal, 'plos2')
			width = 6.83;
		elseif strcmp(journal, 'plos1.5')
			width = 4.86;
		elseif strcmp(journal, 'els1')
			width = 3.543;
		elseif strcmp(journal, 'els2')
			width = 7.48031;
		elseif strcmp(journal, 'els1.5')
			width = 5.51181;
		else
			throw(MException('Argument:Error', 'Invalid journal parameter. See help saveplot'));
		end
		%Set fonts correctly
		if strfind(journal, 'els')
			fs = 7;
		elseif strfind(journal, 'plos')
			fs = 8;
		end
		set(ax,'FontSize',fs, 'FontName', 'Arial');
		set(findall(h,'type','text'),'fontSize',fs,'fontName', 'Arial')
		%ratio = 1.61803398875; %Golden ratio
		ratio = 1.5;
		dimensions = dimensions*[width width/ratio];
	end

	%Add a boundary of 4 pts. I'm not sure if this is working right
	pts = 4*1/72;
	box = [pts, pts, dimensions(1)-2*pts, dimensions(2)-2*pts];

	renderer = '-painters';
	if strcmp(format, 'eps')
		dev = '-depsc';
	elseif strcmp(format, 'png')
		dev = '-dpng';
		renderer = '-zbuffer';
	elseif strcmp(format, 'svg')
		dev = '-dsvg';
		%renderer = '-zbuffer';
	elseif strcmp(format, 'pdf')
		dev = '-dpdf';
	elseif strcmp(format, 'jpg')
		dev = '-djpeg';
		%renderer = '-opengl';
	else
		throw(MException('ArgumentError','Format must be one of eps or png or pdf or jpeg'));
	end

	set(h, 'paperunits', 'inches');
	set(h, 'papersize', dimensions);
	set(h, 'paperposition', box);
	set(h, 'paperpositionmode', 'manual');
	print(h, dev, '-r300', renderer, filename);
end
