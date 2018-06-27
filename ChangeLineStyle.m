function ChangeLineStyle(fh)
TempPath = getapplicationdatadir(fullfile('\Temp'),false,true); 
TempFile = fullfile(TempPath,'tempmatlab.mat');
if exist(TempFile,'file')
    delete(TempFile)
end

ID = 1;
LS = 2;
LW = 4;
MS = 3;
MW = 5;
CL = 6;
hToolbar = findall(fh,'tag','FigureToolBar');
LastItemHandle = findall(hToolbar,'tag','Plottools.PlottoolsOn');
%LastItemHandle.Separator = 'on';
% Add undo dropdown list to the toolbar
pause(0.1)
jToolbar = get(get(hToolbar,'JavaContainer'),'ComponentPeer');

if ~isempty(jToolbar)
    L=findall(fh,'type','line');
    if numel(L) == 0, return, end
    for il = 1 :numel(L)
        Line(il).ID = strcat('Line: ',num2str(il));
        Line(il).LineStyle = set(L(il),'LineStyle');
        Line(il).MarkerStyle = set(L(il),'Marker');
        Line(il).LineWidth = 1:30;
        Line(il).MarkerSize = 1:30;
        Line(il).Color = {L(il).Color,'r','g','b','c','m','y','k','w'};
    end
    
    jCombo(ID) = javax.swing.JComboBox({Line.ID});%(undoActions(end:-1:1));
    set(jCombo(ID),'Name','ID')
    jCombo(ID).setMaximumSize(java.awt.Dimension(100,31));
    jCombo(ID).setMinimumSize(java.awt.Dimension(100,31));
    
    jCombo(LS) = javax.swing.JComboBox;
    set(jCombo(LS),'Name','LS')
    jCombo(LS).setMaximumSize(java.awt.Dimension(30,31));
    jCombo(LS).setMinimumSize(java.awt.Dimension(30,31));
    
    jCombo(MS) = javax.swing.JComboBox;%(undoActions(end:-1:1));
    set(jCombo(MS),'Name','MS')
    jCombo(MS).setMaximumSize(java.awt.Dimension(100,31));
    jCombo(MS).setMinimumSize(java.awt.Dimension(100,31));
    
    jCombo(LW) = javax.swing.JComboBox;
    set(jCombo(LW),'Name','LW')
    jCombo(LW).setMaximumSize(java.awt.Dimension(100,31));
    jCombo(LW).setMinimumSize(java.awt.Dimension(100,31));
    
    jCombo(MW) = javax.swing.JComboBox;
    set(jCombo(MW),'Name','MW')
    jCombo(MW).setMaximumSize(java.awt.Dimension(100,31));
    jCombo(MW).setMinimumSize(java.awt.Dimension(100,31));
    
    jCombo(CL) = javax.swing.JComboBox;
    set(jCombo(CL),'Name','CL')
    jCombo(CL).setMaximumSize(java.awt.Dimension(100,31));
    jCombo(CL).setMinimumSize(java.awt.Dimension(100,31));
    
    set(jCombo(ID), 'ActionPerformedCallback', {@SelectLine,Line,jCombo,L});
    set(jCombo(LS), 'ActionPerformedCallback', {@SelectLineStyle,Line,jCombo,L});
    set(jCombo(MS), 'ActionPerformedCallback', {@SelectLineStyle,Line,jCombo,L});
    set(jCombo(LW), 'ActionPerformedCallback', {@SelectLineStyle,Line,jCombo,L});
    set(jCombo(MW), 'ActionPerformedCallback', {@SelectLineStyle,Line,jCombo,L});
    set(jCombo(CL), 'ActionPerformedCallback', {@SelectLineStyle,Line,jCombo,L});
    
    NC=get(jToolbar,'Components');
    jToolbar(1).add(jCombo(ID),numel(NC));
    NC=get(jToolbar,'Components');
    jToolbar(1).add(jCombo(LS),numel(NC));
    NC=get(jToolbar,'Components');
    jToolbar(1).add(jCombo(LW),numel(NC));
    NC=get(jToolbar,'Components');
    jToolbar(1).add(jCombo(CL),numel(NC));
    NC=get(jToolbar,'Components');
    jToolbar(1).add(jCombo(MS),numel(NC));
    NC=get(jToolbar,'Components');
    jToolbar(1).add(jCombo(MW),numel(NC));
    
    jToolbar(1).repaint;
    jToolbar(1).revalidate;

icon = fullfile(matlabroot,'/toolbox/matlab/icons/tool_shape_stroke.png');
[cdata,map] = imread(icon,'png');
cdata = im2uint8(cdata);
cdata(find(cdata == 0)) = NaN;
hToolbar = findall(fh,'tag','FigureToolBar');
hUndo = uipushtool('parent',hToolbar,'cdata',cdata,'tooltip','Save','ClickedCallback',{@SaveConfig,Line,jCombo,L});
hUndo.Separator = 'on';
end



end

function SaveConfig(src,event,Line,jCombo,L)
TempPath = getapplicationdatadir(fullfile('\Temp'),false,true);
TempFile = fullfile(TempPath,'tempmatlab.mat');
for il=1:numel(L)
    strcat()
end
save(TempFile,'L');
assignin('base','CurveHandle',L)
end

function SelectLine(hCombo,~,Line,jCombo,~)
itemIndex = get(hCombo,'SelectedIndex')+1;
ID = 1;
LS = 2;
LW = 4;
MS = 3;
MW = 5;
CL = 6;
TempPath = getapplicationdatadir(fullfile('\Temp'),false,true);
TempFile = fullfile(TempPath,'tempmatlab.mat');
if exist(TempFile,'file')
    load(TempFile)
    return
end
jCombo(LS).removeAllItems;
NumElem=numel(Line(itemIndex).LineStyle);
for in = 1:NumElem
    jCombo(LS).addItem(Line(itemIndex).LineStyle{in})
end
jCombo(MS).removeAllItems;
NumElem=numel(Line(itemIndex).MarkerStyle);
for in = 1:NumElem
    jCombo(MS).addItem(Line(itemIndex).MarkerStyle{in})
end
NumElem=numel(Line(itemIndex).LineWidth);
for in = 1:NumElem
    jCombo(LW).addItem(num2str(Line(itemIndex).LineWidth(in),'%d'))
end
NumElem=numel(Line(itemIndex).MarkerSize);
for in = 1:NumElem
    jCombo(MW).addItem(num2str(Line(itemIndex).MarkerSize(in),'%d'))
end
NumElem=numel(Line(itemIndex).Color);
for in = 1:NumElem
    if isnumeric(Line(itemIndex).Color{in})
        jCombo(CL).addItem(['[',num2str(Line(itemIndex).Color{in}),']'])
    else
        jCombo(CL).addItem(Line(itemIndex).Color{in})
    end
    
end

end

function SelectLineStyle(hCombo,~,~,jCombo,L)
ID = 1;
LS = 2;
LW = 4;
MS = 3;
MW = 5;
CL = 6;
IDLine = get(jCombo(ID),'SelectedIndex')+1;
Info = get(hCombo,'SelectedItem');
switch get(hCombo,'Name')
    case 'LS'
        String2Eval = 'LineStyle';
    case 'MS'
        String2Eval = 'Marker';
    case 'LW'
        String2Eval = 'LineWidth';
        Info = str2double(Info);
    case 'MW'
        String2Eval = 'MarkerSize';
        Info = str2double(Info);
    case 'CL'
        String2Eval = 'Color';
end
set(L(IDLine),String2Eval,Info);

end
