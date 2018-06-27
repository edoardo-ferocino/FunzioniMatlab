function varargout = PlotSolus(varargin)
% PLOTSOLUS MATLAB code for PlotSolus.fig
%      PLOTSOLUS, by itself, creates a new PLOTSOLUS or raises the existing
%      singleton*.
%
%      H = PLOTSOLUS returns the handle to a new PLOTSOLUS or the handle to
%      the existing singleton*.
%
%      PLOTSOLUS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTSOLUS.M with the given input arguments.
%
%      PLOTSOLUS('Property','Value',...) creates a new PLOTSOLUS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlotSolus_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlotSolus_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlotSolus

% Last Modified by GUIDE v2.5 21-Apr-2017 15:48:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PlotSolus_OpeningFcn, ...
    'gui_OutputFcn',  @PlotSolus_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end


if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before PlotSolus is made visible.
function PlotSolus_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlotSolus (see VARARGIN)

% Choose default command line output for PlotSolus

handles.output = hObject;
clc
handles=FindAllHandles(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlotSolus wait for user response (see UIRESUME)
% uiwait(handles.PlotSolus);

% --- Outputs from this function are returned to the command line.
function varargout = PlotSolus_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over LoadData.
function LoadData_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TabHandle=FindTabHandle(handles);


iDV = 1;
iCV = 1;

if exist('tempsolus.mat','file')
    load('tempsolus.mat','PathName')
    [FileName,PathName,FilterIndex] = uigetfile('*.mat','Load Sim.mat',PathName);
    if FilterIndex ==0 , return, end
    ClearTable(TabHandle,handles)
else
    [FileName,PathName,FilterIndex] = uigetfile('*.mat','Load Sim.mat');
    if FilterIndex ==0 , return, end
end
save('tempsolus.mat','PathName')
load([PathName,FileName],'P','F','Root','iv','iP','datetxt','geom','type')
if ~isfield(P,'Help')
    P(1).Help = [];
end
if ~isfield(F,'Help')
    F(1).Help = [];
end
set(hObject,'String',[PathName,FileName])
%%%
% ODP(1).Name = 'geom';
% ODP(1).Value = geom;
% ODP(1).Unit = '';
% ODP(2).Name = 'type';
% ODP(2).Value = type;
% ODP(2).Unit = '';
%%%
% DefaultVar = zeros(iP-iv);
% ChangingVar = zeros(iv);
for iP = 1:numel(P)
    if(P(iP).Order==0||P(iP).Order == -1)
        Variable(handles.Default).PosInP(iDV) = iP;
        DefaultVar(iDV)=iP;
        TabHandle(handles.Default).RowName(iDV) = {['D' num2str(iDV)]};
        TabHandle(handles.Default).Data(iDV,handles.Name) = {P(iP).Label};
        if mod(P(iP).Default,1)==0
                format = '%d';
        else
                format = '%g';
        end
        TabHandle(handles.Default).Data(iDV,handles.Value) = {num2str(P(iP).Default,format)};
        TabHandle(handles.Default).Data(iDV,handles.Unit) = {P(iP).Unit};
        iDV = iDV +1;
    else
        Variable(handles.Changing).PosInP(P(iP).Order) = iP;
        ChangingVar(P(iP).Order)=iP;
        TabHandle(handles.Changing).RowName(P(iP).Order) = {['C' num2str(P(iP).Order)]};
        TabHandle(handles.Changing).Data(P(iP).Order,handles.Name) = {P(iP).Label};
        if mod(P(iP).Range,1)==0
                format = '%d;';
        else
                format = '%g;';
        end
        TabHandle(handles.Changing).Data(P(iP).Order,handles.Value) = {num2str(P(iP).Range,format)};
        TabHandle(handles.Changing).Data(P(iP).Order,handles.Unit) = {P(iP).Unit};
        iCV = iCV +1;
    end
end
%%%%
% for iOD = 1:numel(ODP)
%     DefaultVar(iDV-1+iOD)=numel(P)+1;
%     P=addP(-1,numel(P)+1,ODP(iOD).Value,0,ODP(iOD).Name,ODP(iOD).Unit,0,P);
%     TabHandle(handles.Default).RowName(iDV+iOD-1) = {['D' num2str(iDV+iOD-1)]};
%     TabHandle(handles.Default).Data(iDV+iOD-1,handles.Name) = {ODP(iOD).Name};
%     TabHandle(handles.Default).Data(iDV+iOD-1,handles.Value) = {ODP(iOD).Value};
% end
%%%%
for iF = 1:numel(F)
    TabHandle(handles.FiguresOfMerit).RowName(iF) = {['F' num2str(iF)]};
    TabHandle(handles.FiguresOfMerit).Data(iF,handles.Name) = {F(iF).Label};
    if mod(F(iF).Levels,1)==0
                format = '%d;';
        else
                format = '%g;';
    end
    TabHandle(handles.FiguresOfMerit).Data(iF,handles.Value) = {num2str(F(iF).Levels,format)};
    TabHandle(handles.FiguresOfMerit).Data(iF,handles.Unit) = {F(iF).Unit};
end

NCV = iCV-1; %Number Changing Variables
ColoumnsFoM = {F.Label};
TabHandle(handles.SelectFoMContour).ColumnName = ColoumnsFoM;
TabHandle(handles.SelectFoMPlot).ColumnName = ColoumnsFoM;
ContPath=CreatePath(PathName,'C',ColoumnsFoM,P,ChangingVar);
PlotPath=CreatePath(PathName,'P',ColoumnsFoM,P,ChangingVar);
for iF = 1:numel(F)
    if NCV>=5
        for iR=1:length(P(ChangingVar(NCV)).Range)
            TabHandle(handles.SelectFoMContour).RowName(iR) = {[P(ChangingVar(NCV)).Label '=' num2str(P(ChangingVar(NCV)).Range(iR))]};
            TabHandle(handles.SelectFoMPlot).RowName(iR) = {[P(ChangingVar(NCV)).Label '=' num2str(P(ChangingVar(NCV)).Range(iR))]};
        end
    end
end

TabHandle(handles.SelectFoMContour).Data = ContPath;
TabHandle(handles.SelectFoMPlot).Data = PlotPath;
handles.PathName = PathName;
handles.FileName = FileName;
handles.NCV = NCV;
handles.P = P;
handles.Coloumns = ColoumnsFoM;
handles.ChangingVar = ChangingVar;
handles.DefaultVar = DefaultVar;
handles.Variable = Variable;
handles.F = F;
handles.TabHandle = TabHandle;
guidata(hObject,handles);

function [NewPath] = CreatePath(Path,PlotType,FoM_Name,P,ChangingVar)
NCV = numel(ChangingVar);
Format = get(findobj('Tag','Format'),'Value');
FormatList = get(findobj('Tag','Format'),'String');
Format = FormatList{Format};
BsIndex = strfind(Path,'\');
Root = Path(BsIndex(end-1)+1:BsIndex(end)-1);
Root=[Path,Root];
if ~iscell(FoM_Name)
    FoM_Name = {FoM_Name};
end
if NCV>=5
    NewPath = cell(length(P(ChangingVar(NCV)).Range),numel(FoM_Name));
else
    NewPath = cell(1,numel(FoM_Name));
end

for iFoM = 1:numel(FoM_Name)
    if NCV>=5
        for iNCV = 1:length(P(ChangingVar(NCV)).Range)
            NewPath{iNCV,iFoM} = strcat(Root,'_',num2str(iNCV),'_',PlotType,'_',FoM_Name{iFoM},Format);
        end
    else
        NewPath{1,iFoM} = strcat(Root,'_',PlotType,'_',FoM_Name{iFoM},Format);
    end
end

% --- Executes when selected cell(s) is changed in FoMSelectContour.
function FoMSelectContour_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to FoMSelectContour (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
OpenPath(hObject,eventdata,handles)

% --- Executes when selected cell(s) is changed in FoMSelectPlot.
function FoMSelectPlot_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to FoMSelectPlot (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
OpenPath(hObject,eventdata,handles)

% --- Executes on selection change in Format.
function Format_Callback(hObject, eventdata, handles)
% hObject    handle to Format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Value')
    case 1
        Fileformat = '.pdf';
    case 2
        Fileformat = '.jpg';
    case 3
        Fileformat = '.fig';
end
if isstruct(handles)
    ContPath=CreatePath(handles.PathName,'C',handles.Coloumns,handles.P,handles.ChangingVar);
    PlotPath=CreatePath(handles.PathName,'P',handles.Coloumns,handles.P,handles.ChangingVar);
    handles.TabHandle(handles.SelectFoMContour).Data = ContPath;
    handles.TabHandle(handles.SelectFoMPlot).Data = PlotPath;
    handles.FileFormat = Fileformat;
end
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns Format contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Format

% --- Executes during object creation, after setting all properties.
function Format_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in GenerateReport.
function GenerateReport_Callback(hObject, eventdata, handles)
% hObject    handle to GenerateReport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(exist([handles.PathName 'Report.txt'],'file'))
    [FileName,PathName,FilterIndex] = uiputfile('*.txt','Type Name of the Report',[handles.PathName 'Report']);
    if FilterIndex == 0, return, end
    fid=fopen([PathName FileName],'w+');
else
    FilterIndex=0;
    fid=fopen([handles.PathName 'Report.txt'],'w+');
end
fprintf(fid,'%s\r\n','---Fixed Variables---');
for iFV = 1:length(handles.DefaultVar)
    if(handles.P(handles.DefaultVar(iFV)).Order == -1)%%% TEMP
        fprintf(fid,'%s\t%s\t%s\r\n',handles.P(handles.DefaultVar(iFV)).Label,handles.P(handles.DefaultVar(iFV)).Default,handles.P(handles.DefaultVar(iFV)).Unit);
    else
        if mod(handles.P(handles.DefaultVar(iFV)).Default,1)==0
                %format = '%s\t%d\t%s\r\n';
                format = '%s\t%d\t%s\t%s\r\n';
        else
                %format = '%s\t%g\t%s\r\n';
                format = '%s\t%g\t%s\t%s\r\n';
        end
%         fprintf(fid,format,handles.P(handles.DefaultVar(iFV)).Label,handles.P(handles.DefaultVar(iFV)).Default,handles.P(handles.DefaultVar(iFV)).Unit);
        S=size(handles.P(handles.DefaultVar(iFV)).Help);
        if(S(1)>0)
            buffer =handles.P(handles.DefaultVar(iFV)).Help(1,:);
        for iRow = 2:S(1)
             buffer=char(strcat(buffer,{' '},handles.P(handles.DefaultVar(iFV)).Help(iRow,:)));
        end
        else
            buffer = handles.P(handles.DefaultVar(iFV)).Help;
        end
        fprintf(fid,format,handles.P(handles.DefaultVar(iFV)).Label,handles.P(handles.DefaultVar(iFV)).Default,handles.P(handles.DefaultVar(iFV)).Unit,buffer);
    end
end
fprintf(fid,'\r\n%s\r\n','---Study Variables---');
for iCV = 1:length(handles.ChangingVar)
    fprintf(fid,'%s\t',handles.P(handles.ChangingVar(iCV)).Label);
    if mod(handles.P(handles.ChangingVar(iCV)).Range,1)==0
                format = '%d ';
        else
                format = '%g ';
    end
    fprintf(fid,'%s\t',num2str(handles.P(handles.ChangingVar(iCV)).Range,format));
    S=size(handles.P(handles.ChangingVar(iCV)).Help);
        if(S(1)>0)
            buffer =handles.P(handles.ChangingVar(iCV)).Help(1,:);
        for iRow = 2:S(1)
             buffer=char(strcat(buffer,{' '},handles.P(handles.ChangingVar(iCV)).Help(iRow,:)));
        end
        else
            buffer = handles.P(handles.ChangingVar(iCV)).Help;
        end
    fprintf(fid,'%s\t',handles.P(handles.ChangingVar(iCV)).Unit);
    fprintf(fid,'%s\r\n',buffer);
end
fprintf(fid,'\r\n%s\r\n','---Figures of Merit---');
for iFoM = 1:length(handles.F)
    fprintf(fid,'%s\t',handles.F(iFoM).Label);
    if mod(handles.F(iFoM).Levels,1)==0
                format = '%d ';
        else
                format = '%g ';
    end
    fprintf(fid,'%s\t',num2str(handles.F(iFoM).Levels,format));
    fprintf(fid,'%s\t',handles.F(iFoM).Unit);
    S=size(handles.F(iFoM).Help);
        if(S(1)>0)
            buffer = handles.F(iFoM).Help(1,:);
        for iRow = 2:S(1)
             buffer=char(strcat(buffer,{' '},handles.F(iFoM).Help(iRow,:)));
        end
        else
            buffer = handles.F(iFoM).Help;
        end
    fprintf(fid,'%s\r\n',buffer);
end
if isfield(handles,'NoteText')
    fprintf(fid,'\r\n---Note---\r\n');
    fprintf(fid,'%s\r\n',string(handles.NoteText));
end
fid=fclose(fid);
if fid==0
    if FilterIndex==0
        msgbox({'Report created at:',handles.PathName,'as:','Report.txt'},'Done!','modal');
    else
        msgbox({'Report created at:',handles.PathName,'as:',FileName},'Done!','modal');
    end
    
else
    msgbox('Report NOT created','ERROR!','modal');
end

% --- Executes when user attempts to close PlotSolus.
function PlotSolus_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to PlotSolus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if exist('tempsolus.mat','file')
    delete('tempsolus.mat')
end
% Hint: delete(hObject) closes the figure
delete(hObject);

function Note_Callback(hObject, eventdata, handles)
% hObject    handle to Note (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.NoteText = char(hObject.String);
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of Note as text
%        str2double(get(hObject,'String')) returns contents of Note as a double

% --- Executes during object creation, after setting all properties.
function Note_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Note (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected cell(s) is changed in DefaultParameters.
function DefaultParameters_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to DefaultParameters (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
FillDescriptor(hObject,eventdata,handles)

% --- Executes when selected cell(s) is changed in ChangingParameters.
function ChangingParameters_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to ChangingParameters (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
FillDescriptor(hObject,eventdata,handles)

% --- Executes when selected cell(s) is changed in FiguresOfMerit.
function FiguresOfMerit_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to FiguresOfMerit (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
FillDescriptor(hObject,eventdata,handles)

function FillDescriptor(hObject,eventdata,handles)
if ~isempty(eventdata.Indices)
if isfield(handles.P,'Help')
    if strcmpi(hObject.Tag,'DefaultParameters')
        HelpString = handles.P(handles.DefaultVar(eventdata.Indices(1))).Help;
        IS_P_HELP_FIELD = 1;
    else
        if strcmpi(hObject.Tag,'ChangingParameters')
            HelpString = handles.P(handles.ChangingVar(eventdata.Indices(1))).Help;
            IS_P_HELP_FIELD = 1;
        else
            IS_P_HELP_FIELD = 0;
        end
    end
else
    IS_P_HELP_FIELD = 0;
end
if isfield(handles.F,'Help')
    if strcmpi(hObject.Tag,'FiguresOfMerit')
        HelpString = handles.F(eventdata.Indices(1)).Help;
        IS_F_HELP_FIELD = 1;
    else
        IS_F_HELP_FIELD = 0;
    end
else
    IS_F_HELP_FIELD = 0;
end
IS_HELP_FIELD = IS_P_HELP_FIELD || IS_F_HELP_FIELD;
if IS_HELP_FIELD && ~isempty(HelpString)
    set(findobj('Tag','Descriptor'),'String',HelpString)
else
    set(findobj('Tag','Descriptor'),'String','Help')
end
handles.IS_P_HELP_FIELD = IS_P_HELP_FIELD;
handles.IS_F_HELP_FIELD = IS_F_HELP_FIELD;
guidata(hObject,handles)
end

function OpenPath(hObject,eventdata,handles)
if ~isempty(eventdata.Indices)
    Path = hObject.Data(eventdata.Indices(1),eventdata.Indices(2));
    try
        msgHandle = msgbox('Opening...');
        winopen(char(Path));
        pause(1);
        delete(msgHandle)
    catch ME
        msgbox(ME.message,'Error')
        uiwait()
    end
end

function ClearTable(TabHandles,handles)
nTB = numel(TabHandles);
EmptyVal = {''};
for it=1:nTB
    TabHandles(it).Data = EmptyVal;
    TabHandles(it).RowName = EmptyVal;
    if it == handles.SelectFoMContour || it== handles.SelectFoMPlot
        TabHandles(it).ColumnName = EmptyVal;
    end
end
NoteHandle=findobj('Tag','Note');
NoteHandle.String='Insert Notes';
DescriptorHandle=findobj('Tag','Descriptor');
DescriptorHandle.String='Help';


function P = addP(Order,iP,Default,Range,Label,Unit,Title,P)

% addP(Order,iP,Default,Range,Label,Unit,Title,P): Add one entry in the parameter space
%
%   Order = Order in the Output (0=take default, 1=x-axis, 2=yaxis, 3=rows, 4=columns
%   iP = index of the P entry
%   Default = default value
%   Range = range of values
%   Label = label used for rapresentation
%   Unit = measurement unit of the value
%   Dim = number of elements = 1 if Order = 0, else num elements in Range
%   P = P structure

P(iP).Order=Order;
if(isempty(P(iP).Default)==0), disp('Error in addP: Duplicated iP'); end
P(iP).Default=Default;
P(iP).Range=Range;
P(iP).Label=Label;
P(iP).Unit=Unit;
P(iP).Title=Title;
nP=numel(P);
if(nP~=iP), disp('Error in addP: Mismatch beetween nP and iP'); end
for i=1:nP
    if(isempty(P(i).Default)==1), disp('Error in addP: Empty Elements'); end
end
if Order==0, P(iP).Dim=1; else P(iP).Dim=numel(Range); end

% --- Executes on button press in LoadReport.
function LoadReport_Callback(hObject, eventdata, handles)
% hObject    handle to LoadReport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TabHandle=FindTabHandle(handles);

if exist('tempsolus.mat','file')
    load('tempsolus.mat','PathName')
    [FileName,PathName,FilterIndex] = uigetfile('*.txt','Load Report.mat',PathName);
    if FilterIndex ==0 , return, end
    ClearTable(TabHandle,handles)
else
    [FileName,PathName,FilterIndex] = uigetfile('*.txt','Load Report.mat');
    if FilterIndex ==0 , return, end
end
save('tempsolus.mat','PathName')
fid=fopen([PathName FileName],'r');
iVar=0;
VarType = 1;
buffer=fgetl(fid);
while(~feof(fid))
buffer=fgetl(fid);
if (~isempty(buffer))
    if buffer(1:3)~='---'
        iVar = iVar +1;
    end
else
  VarTypeNum(VarType) = iVar;
  iVar = 0;
  VarType = VarType +1;
end
end
VarTypeNum(VarType) = iVar;
frewind(fid);
fgetl(fid);
%A=textscan(fid,'%s %s %s',VarTypeNum(handles.Default),'delimiter','\t');
A=textscan(fid,'%s %s %s %s',VarTypeNum(handles.Default),'delimiter','\t');
LabelA=A{1};
ValueA=A{2};
UnitA=A{3};
HelpA=A{4};
for iDV = 1:numel(LabelA)
        DefaultVar(iDV)=iDV;
        TabHandle(handles.Default).RowName(iDV) = {['D' num2str(iDV)]};
        P(iDV).Label = LabelA{iDV};
        TabHandle(handles.Default).Data(iDV,handles.Name) = {P(iDV).Label};
        TabHandle(handles.Default).Data(iDV,handles.Value) = {ValueA{iDV}};
        if isnan(str2double(ValueA{iDV}))
            P(iDV).Default = ValueA{iDV};
            P(iDV).Order = -1;
        else
            P(iDV).Default = str2double(ValueA{iDV});
            P(iDV).Order = 0;
        end
        P(iDV).Unit = UnitA{iDV};
        P(iDV).Help = HelpA{iDV};
        TabHandle(handles.Default).Data(iDV,handles.Unit) = {P(iDV).Unit};
       
end
buffer = fgetl(fid);
buffer = fgetl(fid);
buffer = fgetl(fid);
%A=textscan(fid,'%s %s %s',VarTypeNum(handles.Changing),'delimiter','\t');
A=textscan(fid,'%s %s %s %s',VarTypeNum(handles.Changing),'delimiter','\t');
LabelA=A{1};
RangeA=A{2};
UnitA=A{3};
HelpA=A{4};
for iCV = 1:numel(LabelA)
        ChangingVar(iCV)=iDV+iCV;
        TabHandle(handles.Changing).RowName(iCV) = {['C' num2str(iCV)]};
        P(iCV+iDV).Label = LabelA{iCV};
        TabHandle(handles.Changing).Data(iCV,handles.Name) = {P(iCV+iDV).Label};
        TabHandle(handles.Changing).Data(iCV,handles.Value) = {RangeA{iCV}};%{num2str(P(iCV+iDV).Range)};
        P(iCV+iDV).Range = str2num(RangeA{iCV});
        P(iCV+iDV).Unit = UnitA{iCV};
        TabHandle(handles.Changing).Data(iCV,handles.Unit) = {P(iCV+iDV).Unit};
        P(iCV+iDV).Order = iCV;
        P(iCV+iDV).Help = HelpA{iCV};
end
buffer = fgetl(fid);
buffer = fgetl(fid);
buffer = fgetl(fid);
%A=textscan(fid,'%s %s %s',VarTypeNum(handles.FiguresOfMerit),'delimiter','\t');
A=textscan(fid,'%s %s %s %s',VarTypeNum(handles.FiguresOfMerit),'delimiter','\t');
LabelA=A{1};
RangeA=A{2};
UnitA=A{3};
HelpA=A{4};
for iFoM = 1:numel(LabelA)
        TabHandle(handles.FiguresOfMerit).RowName(iFoM) = {['F' num2str(iFoM)]};
        F(iFoM).Label = LabelA{iFoM};
        TabHandle(handles.FiguresOfMerit).Data(iFoM,handles.Name) = {F(iFoM).Label};
        TabHandle(handles.FiguresOfMerit).Data(iFoM,handles.Value) = {RangeA{iFoM}};
        F(iFoM).Levels = str2num(RangeA{iFoM});
        F(iFoM).Unit = UnitA{iFoM};
        F(iFoM).Help = HelpA{iFoM};
        TabHandle(handles.FiguresOfMerit).Data(iFoM,handles.Unit) = {F(iFoM).Unit};
end
buffer = fgetl(fid);
buffer = fgetl(fid);
buffer = fgetl(fid);
ib=1;
buffer2 = {'Insert Notes'};
while(~feof(fid))
    buffer = fgetl(fid);
    buffer2{ib,1} = buffer;
    ib=ib+1;
end
handles.Note.String = char(buffer2);
Coloumns={F.Label};
NCV=length(ChangingVar);
ContPath=CreatePath(PathName,'C',Coloumns,P,ChangingVar);
PlotPath=CreatePath(PathName,'P',Coloumns,P,ChangingVar);
for iF = 1:numel(F)
    if NCV>=5
        for iR=1:length(P(ChangingVar(NCV)).Range)
            TabHandle(handles.SelectFoMContour).RowName(iR) = {[P(ChangingVar(NCV)).Label '=' num2str(P(ChangingVar(NCV)).Range(iR))]};
            TabHandle(handles.SelectFoMPlot).RowName(iR) = {[P(ChangingVar(NCV)).Label '=' num2str(P(ChangingVar(NCV)).Range(iR))]};
        end
    end
end
TabHandle(handles.SelectFoMContour).ColumnName = Coloumns;
TabHandle(handles.SelectFoMPlot).ColumnName = Coloumns;
TabHandle(handles.SelectFoMContour).Data = ContPath;
TabHandle(handles.SelectFoMPlot).Data = PlotPath;
handles.PathName = PathName;
handles.NCV = NCV;
handles.P = P;
handles.Coloumns = Coloumns;
handles.ChangingVar = ChangingVar;
handles.DefaultVar = DefaultVar;
handles.F = F;
handles.TabHandle = TabHandle;
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of LoadReport

function TabHandle = FindTabHandle(handles)
TabHandle(handles.Default) = findobj('Tag','DefaultParameters');
TabHandle(handles.Changing) = findobj('Tag','ChangingParameters');
TabHandle(handles.FiguresOfMerit) = findobj('Tag','FiguresOfMerit');
TabHandle(handles.SelectFoMContour) = findobj('Tag','FoMSelectContour');
TabHandle(handles.SelectFoMPlot) = findobj('Tag','FoMSelectPlot');





function InsertHelp_Callback(hObject, eventdata, handles)
% hObject    handle to InsertHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.InsertHelp = char(hObject.String);
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of InsertHelp as text
%        str2double(get(hObject,'String')) returns contents of InsertHelp as a double


% --- Executes during object creation, after setting all properties.
function InsertHelp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InsertHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Descriptor.
function Descriptor_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Descriptor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hObject.Visible = 'off';
InsertHelpHandle = findobj('Tag','InsertHelp');
SaveHelpHandle = findobj('Tag','SaveHelp');
CellAdressHandle = findobj('Tag','CellAdress');
SaveToHandle = findobj('Tag','SaveTo');
CancelHelpHandle = findobj('Tag','CancelHelp');

InsertHelpHandle.Visible = 'on';
SaveHelpHandle.Visible = 'on';
CellAdressHandle.Visible = 'on';
SaveToHandle.Visible = 'on';
CancelHelpHandle.Visible = 'on';
InsertHelpHandle.Position = hObject.Position;


% --- Executes on button press in SaveHelp.
function SaveHelp_Callback(hObject, eventdata, handles)
% hObject    handle to SaveHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InsertHelpHandle = findobj('Tag','InsertHelp');
SaveHelpHandle = findobj('Tag','SaveHelp');
CellAdressHandle = findobj('Tag','CellAdress');
SaveToHandle = findobj('Tag','SaveTo');
DescriptorHandle = findobj('Tag','Descriptor');
CancelHelpHandle = findobj('Tag','CancelHelp');

DescriptorHandle.Visible = 'on';
InsertHelpHandle.Visible = 'off';
SaveHelpHandle.Visible = 'off';
CellAdressHandle.Visible = 'off';
SaveToHandle.Visible = 'off';
CancelHelpHandle.Visible = 'off';
P=handles.P;
F=handles.F;
save([handles.PathName handles.FileName],'P','F','-append');
msgbox({'Saved inserted help in P struct of file:',[handles.PathName handles.FileName]},'Done!','modal');

% Hint: get(hObject,'Value') returns toggle state of SaveHelp



function CellAdress_Callback(hObject, eventdata, handles)
% hObject    handle to CellAdress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Adress = hObject.String;
A=sscanf(char(Adress)   ,'%c%d');
switch char(A(1))
    case 'C'
        handles.P(handles.ChangingVar(A(2))).Help = handles.InsertHelp;
    case 'D'
        handles.P(handles.DefaultVar(A(2))).Help = handles.InsertHelp;
    case 'F'
        handles.F(A(2)).Help = handles.InsertHelp;
    otherwise 
        msgbox('Cell not found correctly','Error')
end
guidata(hObject,handles),
% Hints: get(hObject,'String') returns contents of CellAdress as text
%        str2double(get(hObject,'String')) returns contents of CellAdress as a double


% --- Executes during object creation, after setting all properties.
function CellAdress_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellAdress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CancelHelp.
function CancelHelp_Callback(hObject, eventdata, handles)
% hObject    handle to CancelHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InsertHelpHandle = findobj('Tag','InsertHelp');
SaveHelpHandle = findobj('Tag','SaveHelp');
CellAdressHandle = findobj('Tag','CellAdress');
SaveToHandle = findobj('Tag','SaveTo');
CancelHelpHandle = findobj('Tag','CancelHelp');
DescriptorHandle = findobj('Tag','Descriptor');

DescriptorHandle.Visible = 'on';
InsertHelpHandle.Visible = 'off';
SaveHelpHandle.Visible = 'off';
CellAdressHandle.Visible = 'off';
SaveToHandle.Visible = 'off';
CancelHelpHandle.Visible = 'off';
% Hint: get(hObject,'Value') returns toggle state of CancelHelp

function [handles]=FindAllHandles(handles)
handles.Default = 1;
handles.Changing = 2;
handles.FiguresOfMerit = 3;
handles.SelectFoMContour = 4;
handles.SelectFoMPlot = 5;
handles.Name = 1;
handles.Value = 2;
handles.Unit = 3;
handles.Help = 4;
handles.DefaultLabel = 'D';
handles.ChangingLabel = 'C';
handles.FiguresOfMeritLabel = 'F';
