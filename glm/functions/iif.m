function result = iif(varargin) 
	%iif		Returns values depending on conditionals provided
	%
	% Usage:
	%               result = iif(varargin)
	%
	% Input:
	%               varargin = alternating list of conditionals and corresponding result. See example
	%
	% Examples:
	%               normalize = @(x) iif( ~all(isfinite(x)), @() error('Must be finite!'), ...
        %                 all(x == 0),       @() zeros(size(x)), ...
        %                 true,              @() x/norm(x) );
	%		normalize([1 1 0])
	%		normalize([0 inf 2])

	result = varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
end
