function plotmult(h, filename, nframes, format, dimensions, journal)
%plotmult       Saves the specified plot consisting of nframes of size specified in dimensions
%
% Usage:
%                       plotmult(h, filename, nframes, format, dimensions)
%
% Input:
%                       h = figure object. Can just use 'gcf'
%                       filename = output filename
%                       nframes = number of frames in image
%                       format = one of 'eps', 'png' or 'pdf' (optional, default 'eps')
%                       dimensions = (optional, default = [3.27, 2.18], PloS one column) width and height of figure in inches OR
%                                       a scale factor to multiply the journal specific dims described below
%                       journal = (optional) if specified as one of 'plos1', 'plos2', 'plos1.5' or 'els1', 'els2', 'els1.5'
%                                       will set width to that journals specified one column, two column, 1.5 column width (inches)
%                                       which can then be multiplied by the scale factor in dimensions if this plot is part
%                                       of a subplot -- otherwise just set dimensions to 1. This is a bit redundant for this function
%                                       as to do this properly would require using saveplot...
%
% Examples:
%                       subplot(1,2,1)
%                       plot(0:10, (0:10).^2);
%                       subplot(1,2,2)
%                       plot(0:10, exp(0:10));
%                       plotmult(gcf, './plots/test.eps', 2, 'eps', 1, 'plos1');

        if (nargin < 4)
                format = 'eps';
                dimensions = [3.27 2.18];
        elseif nargin < 5
                dimensions = [3.27 2.18];
        elseif nargin == 5
                if length(dimensions) == 1
                        throw(MException('Argument:Error', 'If specifying a scale dimension, must also specify journal. See help saveplot'));
                end
        elseif nargin == 6
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
                        fs = 10;
                end
                set(ax,'FontSize',fs, 'FontName', 'Arial');
                set(findall(h,'type','text'),'fontSize',fs,'fontName', 'Arial')
                %ratio = 1.61803398875; %Golden ratio
                ratio = 1.5;
                dimensions = dimensions*[width width/ratio];
        end

        dimensions(1) = dimensions(1)*nframes;
        box = [0 0 dimensions];

        if strcmp(format, 'eps')
                dev = '-depsc';
        elseif strcmp(format, 'png')
                dev = '-dpng';
        elseif strcmp(format, 'pdf')
                dev = '-dpdf';
        else
                throw(MException('ArgumentError','Format must be one of eps or png or pdf'));
        end

        set(h, 'paperunits', 'inches');
        set(h, 'papersize', dimensions);
        set(h, 'paperposition', box);
        print(h, filename, dev, '-painters');
end

