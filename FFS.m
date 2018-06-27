function [ varargout ] = FFS( varargin )
W='MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame';
isFH = false;
warning('off',W);
for iv=1:numel(varargin)
    if strcmpi(varargin{iv},'visible') && strcmpi(varargin{iv+1},'off')
        varargin{iv+1} = 'on';
        S=warning('query','backtrace');
        warning('off','backtrace');
        warning('FFS:VisiblePropOff','Visible property must be ''on''');
        warning(S);
    end
    if ishandle(varargin{iv})
        isFH = true;
        FH = varargin{iv};
    end
end
if ~isFH, FH=figure(varargin{:}); end
jFrame = get(handle(FH),'JavaFrame');

while(~jFrame.isMaximized)
    pause(0.5)
    jFrame.setMaximized(true);
    pause(0.5)
end
warning('on',W);
if nargout 
    varargout{1} = FH;
end
end

