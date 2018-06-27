function DefineBorder(fh,PaceY,PaceX)
icon = fullfile(matlabroot,'/toolbox/matlab/icons/plotpicker-artimeseries.png');
[cdata,~] = imread(icon,'png');
hToolbar = findall(fh,'tag','FigureToolBar');
hUndo = uipushtool('parent',hToolbar,'cdata',cdata, 'tooltip','Pick X Points','ClickedCallback',{@AddPickXY_CB,fh,PaceY,PaceX});
hUndo.Separator = 'on';

    function POS=AddPickXY_CB(~,~,fh,PaceY,PaceX)
        POS = [0 0];
        [PX,PY] = getpts(fh);
        PX = round(PX);
        PY = round(PY);
        P = ([PX PY]);
        for iy = 1:numel(PY)-1;
            if PY(iy+1)>= PY(iy);
                PaY = abs(PaceY);
            else
                PaY = -abs(PaceY);
            end
            Yq = PY(iy):PaY:PY(iy+1);
            Yq = Yq';
            Xvint = interp1([PY(iy) PY(iy+1)],[PX(iy) PX(iy+1)],Yq);
            Xvint = round(Xvint);
            POS = [POS; [Xvint Yq]];
        end
        POS = POS(2:end,:);
        assignin('base', 'POS', POS)
        nfh = copyobj(fh,fh.Parent);
        for iy = 1:length(POS)-1
            if POS(iy+1,2)~=POS(iy,2)
                if  POS(iy+1,2)>=POS(iy,2)
                    nfh.Children.Children.CData(POS(iy,2),1:POS(iy,1)) = 0;
                else
                    nfh.Children.Children.CData(POS(iy,2),POS(iy,1):end) = 0;
                end
            end
        end
        Top = max(POS(:,2));
        nfh.Children.Children.CData(Top:end,:) = 0;
        DefineBorder(nfh,PaceY,PaceX);
    end

end