function PickCurve(fh,Data)
icon = fullfile(matlabroot,'/toolbox/matlab/icons/plotpicker-glyphplot.png');
[cdata,~] = imread(icon,'png');
hToolbar = findall(fh,'tag','FigureToolBar');
hUndo = uipushtool('parent',hToolbar,'cdata',cdata, 'tooltip','Pick X Points','ClickedCallback',{@PickCurve_CB,fh});
hUndo.Separator = 'on';

    function PickCurve_CB(~,~,fh)
        if ~isempty(findobj(fh,'type','UIContextMenu'))
            c = findobj(fh,'type','UIContextMenu');
            delete(c)
            msgbox('Function Disabled','AddPickXY OFF')
            return
        else
%             msgh = msgbox('Function Enabled','AddPickXY ON');
%             pause(0.5)
%             delete(msgh);
        end
        figure(fh)
        [x,y] = ginput(1);
        roi = [round(x) round(y)];
        figure(100)
        semilogy(1:numel(Data(roi(2),roi(1),:)),squeeze(Data(roi(2),roi(1),:)))
    end
end
