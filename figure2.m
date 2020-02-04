% [EDITED, 2018-06-05, typos fixed]
function FigHandle = figure2(varargin)
[figureargs,otherargs]=parseInputs(varargin);
isuifigure = otherargs.uifigure;
if ~isempty(figureargs)
    if any(ishandle(figureargs{1}))||any(isnumeric(figureargs{1})), isuifigure = 0; end
end
MP = get(0, 'MonitorPositions');
if size(MP, 1) == 1  % Single monitor
    if isuifigure
        FigH = uifigure(figureargs{:});
    else
        FigH = figure(figureargs{:});
    end
else                 % Multiple monitors
    % Catch creation of figure with disabled visibility:
    indexVisible = find(strncmpi(figureargs(1:2:end), 'Vis', 3));
    if ~isempty(indexVisible)
        paramVisible = figureargs(indexVisible(end) + 1);
    else
        paramVisible = get(0, 'DefaultFigureVisible');
    end
    %
    Shift = MP(2, 1:2);
    try
        if isuifigure
            FigH = uifigure(figureargs{:});
        else
            FigH = figure(figureargs{:});
        end
    catch ME
        throwAsCaller(ME)
    end
    if ~any(ishandle(figureargs{1}))&&~any(isnumeric(figureargs{1}))
        set(FigH, 'Units', 'pixels');
        pos      = get(FigH, 'Position');
        set(FigH, 'Position', [pos(1:2) + Shift, pos(3:4)], ...
            'Visible', paramVisible);
    end
end
if nargout ~= 0
    FigHandle = FigH;
end

end

function [figureargs,otherargsstruc]=parseInputs(argin)
figureargs = cell.empty;
fields={'uifigure'};
val = {false};
for ifl = 1:numel(fields)
    otherargsstruc.(fields{ifl}) = val{ifl};
end
if isempty(argin)
    return;
elseif ishandle(argin{1})
    figureargs{numel(figureargs)+1} = argin{1};
    argin = argin(2:end);
elseif isnumeric(argin{1})
    figureargs{numel(figureargs)+1} = argin{1};
    argin = argin(2:end);
end
tfh={'Position' 'OuterPosition' 'InnerPosition' 'Units' 'Renderer' 'RendererMode' 'Visible' 'Color' 'CloseRequestFcn' 'CurrentAxes' 'CurrentCharacter' 'CurrentObject' 'CurrentPoint' 'SelectionType' 'SizeChangedFcn' 'FileName' 'IntegerHandle' 'NextPlot' 'Alphamap' 'Colormap' 'GraphicsSmoothing' 'WindowButtonDownFcn' 'WindowButtonUpFcn' 'WindowButtonMotionFcn' 'WindowScrollWheelFcn' 'WindowKeyPressFcn' 'WindowKeyReleaseFcn' 'MenuBar' 'ToolBar' 'Pointer' 'PointerShapeCData' 'PointerShapeHotSpot' 'Name' 'NumberTitle' 'Number' 'Children' 'Parent' 'HandleVisibility' 'UIContextMenu' 'ButtonDownFcn' 'BusyAction' 'BeingDeleted' 'Interruptible' 'CreateFcn' 'DeleteFcn' 'Type' 'Tag' 'UserData' 'Clipping' 'Scrollable' 'WindowStyle' 'WindowState' 'DockControls' 'Resize' 'PaperPosition' 'PaperPositionMode' 'PaperSize' 'PaperType' 'PaperUnits' 'InvertHardcopy' 'PaperOrientation' 'KeyPressFcn' 'KeyReleaseFcn'};
isprops = arrayfun(@(genprop) strcmpi(argin,tfh(genprop)),1:numel(tfh),'UniformOutput',false)';
isprops = any(cell2mat(isprops));
isprops(find(isprops == 1)+1)=1;
nfigargs = numel(argin(isprops));
nactfigargs = numel(figureargs);
figureargs(nactfigargs+1:nactfigargs+nfigargs) = argin(isprops);
otherargs = argin(~isprops);
for iot = 1:2:numel(otherargs)
    switch lower(otherargs{iot})
        case 'uifigure'
            val = str2logic(otherargs{iot+1});
    end
    otherargsstruc.(lower(otherargs{iot})) = val;
end
end
function str = str2logic(str)
if strcmpi(str,'true'),str = true;else, str = false;end
end