function h = figureName(varargin)
fh=figure(varargin{:});
if any(strcmpi(varargin,'Name'))
    fh.NumberTitle = 'off';
end
if nargout
    h = fh;
end
end

