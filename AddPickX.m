function AddPickX(fh)
icon = fullfile(matlabroot,'/toolbox/matlab/icons/plotpicker-glyphplot.png');
[cdata,~] = imread(icon,'png');
hToolbar = findall(fh,'tag','FigureToolBar');
hUndo = uipushtool('parent',hToolbar,'cdata',cdata, 'tooltip','Pick X Points','ClickedCallback',{@AddPickXY_CB,fh});
hUndo.Separator = 'on';

    function AddPickXY_CB(~,~,fh)
        if ~isempty(findobj(fh,'type','UIContextMenu'))
            c = findobj(fh,'type','UIContextMenu');
            delete(c)
            msgbox('Function Disabled','AddPickXY OFF')
            return
        else
            %msgbox('Function Enabled','AddPickXY ON')
        end
        c = uicontextmenu(fh);
        lh = findall(fh,'type','line');
        for ic = 1:length(lh)
            lh(ic).UIContextMenu = c;
        end
        NZ = 7;
        for iz = 1:NZ
            uimenu('Parent',c,'Label',strcat('Zone: ',num2str(iz)),'Callback',{@SetZoneCB,iz});
        end
        
        function SetZoneCB(src,~,Zn)
            NM = numel(src.Parent.Children);
            for in = 1:NM
                src.Parent.Children(in).Checked = 'off';
            end
            src.Checked = 'on';
            % hAllAxes = findobj(src.Parent.Parent,'type','axes');
            % hLeg = findobj(hAllAxes,'tag','legend');
            % hAxes = setdiff(hAllAxes,hLeg); % All axes which are not
            [PX,PY]=getpts(src.Parent.Parent);
            P = [PX PY];
            try
                Zone = evalin('base','Zone');
                Zone(Zn).PickPoints(end+1) = {P};
            catch
                Zone(Zn).PickPoints(1) = {P};
            end
            assignin('base', 'Zone', Zone)
        end
        
    end

end
