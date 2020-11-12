function [ Data, varargout ] = DatRead(FileName,varargin)
%DatRead('FileName')
%Can be as input a selection of the following parameters
%DatRead(...,'nSource',nSource,'nDet',nDet,'nBoard',nBoard','loop4',numloop4,'loop5',numloop5,'datatype','uint32'/'ushort'/'double','compilesub',true/false,'forcereading',true/false)
%[Data,Header,EasyReadableHead,SubHeaders,EasyReadableSubHead,UnSqueezedHeader]=DatRead3(...)

%% Get #inputs/outputs
NumArgOut = nargout-1;
NumArgin = nargin-1;

%% Set default values
ForceReading = false;
isCompileSubHeader = false;
HeadLen = 764;
datatype = 'ushort';
NumBoard = 1; NumDet= 1; NumSource = 1;
NumMandArgs = 8;
loop5 = 1;loop4 = 1;
ismandatoryarg = zeros(NumMandArgs,1);
isdatatypeargin=0;
LOOP5 = 1;LOOP4 = 2;LOOP3 = 3;LOOP2 = 4;LOOP1 = 5;SOURCE = 6;DET = 7;BOARD = 8;BIN = 9;
Data = double.empty;
varargout = cell.empty(NumArgOut,0);

%% Check if file exists
if isfile(FileName)
    FilePath = FileName;
else
    FilePath = strcat(FileName,'.DAT');
end
if ~isfile(FilePath)
    errordlg('File not found');
    return;
end

%% Checking inputs parameters
for iN = 1:NumArgin
    switch lower(varargin{iN})
        case 'loop5'
            loop5 = varargin{iN+1};
            ismandatoryarg(LOOP5)=1;
        case 'loop4'
            loop4 = varargin{iN+1};
            ismandatoryarg(LOOP4)=1;
        case 'loop3'
            loop3 = varargin{iN+1};
            ismandatoryarg(LOOP3)=1;
        case 'loop2'
            loop2 = varargin{iN+1};
            ismandatoryarg(LOOP2)=1;
        case 'loop1'
            loop1 = varargin{iN+1};
            ismandatoryarg(LOOP1)=1;
        case 'datatype'
            datatype = varargin{iN+1};
            switch datatype
                case 'ushort'
                    datasize = 2;
                case 'uint32'
                    datasize = 4;
                case 'double'
                    datasize = 8;
            end
            isdatatypeargin = 1;
        case 'compilesub'
            isCompileSubHeader = varargin{iN+1};
        case 'forcereading'
            if(ischar(varargin{iN+1})||isstring(varargin{iN+1}))
                varargin{iN+1}=string2boolean(varargin{iN+1});
            end
            ForceReading = logical(varargin{iN+1});
        case 'nsource'
            NumSource = varargin{iN+1};
            ismandatoryarg(SOURCE)=1;
        case 'ndet'
            NumDet = varargin{iN+1};
            ismandatoryarg(DET)=1;
        case 'nboard'
            NumBoard = varargin{iN+1};
            ismandatoryarg(BOARD)=1;
        case 'nbin'
            NumBin = varargin{iN+1};
            ismandatoryarg(BIN)=1;
    end
end

%% Check if all the necessary inputs are present
if any(ismandatoryarg([BOARD DET SOURCE])) == true ...
        && sum(ismandatoryarg([BOARD DET SOURCE]))~=numel(ismandatoryarg([BOARD DET SOURCE])) == true
    errordlg('Not all the necessary inputs. Please insert nDet, nSource, nBoard or none of them');
    return;
end

fid=fopen(FilePath,'rb');
Head=fread(fid,HeadLen,'uint8');
if numel(unique(Head)) ==1 ...
        && sum(ismandatoryarg)~=numel(ismandatoryarg) == true
    errordlg('Please insert all loop values and nDet, nSource, nBoard, NumBin');
    fclose(fid);
    return;
end

CompiledHeader = FillHeader(Head);
SubLen=CompiledHeader.SizeSubHeader;

if CompiledHeader.SubHeader == 0, SubLen = 0;end
if SubLen == 0 && numel(unique(Head))==1, SubLen = 204; end
if SubLen == 0
    SkipSub = true;
    if(~all(ismandatoryarg([BOARD DET SOURCE])))
        errordlg('Subheader missing: Please insert nSource, nDet, nBoard');
        fclose(fid);
        return;
    end
else
    SkipSub = false;
end

%% Assign updated values
if ~exist('NumBin','var'), NumBin = CompiledHeader.McaChannNum;end
if ~exist('datasize','var'), datasize = CompiledHeader.SizeData; end
CompiledHeader.LoopNum(4) = loop4; CompiledHeader.LoopNum(5) = loop5;
CompiledHeader.LoopFirst(4) = 0; CompiledHeader.LoopFirst(5) = 0;
CompiledHeader.LoopLast(4) = 0; CompiledHeader.LoopLast(5) = 0;
CompiledHeader.LoopDelta(4) = 1; CompiledHeader.LoopDelta(5) = 1;
CompiledHeader.LoopHome(4) = 0; CompiledHeader.LoopHome(5) = 0;
if exist('loop3','var'), CompiledHeader.LoopNum(3) = loop3; end
if exist('loop2','var'), CompiledHeader.LoopNum(2) = loop2; end
if exist('loop1','var'), CompiledHeader.LoopNum(1) = loop1; end
CompiledHeader.LoopNum(1:3)=flip(CompiledHeader.LoopNum(1:3));
CompiledHeader.LoopFirst(1:3)=flip(CompiledHeader.LoopFirst(1:3));
CompiledHeader.LoopLast(1:3)=flip(CompiledHeader.LoopLast(1:3));
CompiledHeader.LoopDelta(1:3)=flip(CompiledHeader.LoopDelta(1:3));
CompiledHeader.LoopHome(1:3)=flip(CompiledHeader.LoopHome(1:3));

%% Get number of sources,detector and boards
DataSize = [2,4,8];
DataType = {'ushort','long','double'};
if isdatatypeargin == true
    DataType = {datatype};
    DataSize = datasize;
end
info=dir(FilePath);
for itry = 1:numel(DataType)
    if SkipSub == false
        [NumSource,NumDet,NumBoard]=GetSourceDetBoard(fid,HeadLen,SubLen,NumBin,DataType{itry});
    end
    NumLoop=CompiledHeader.LoopNum;
    if info.bytes == (HeadLen + prod(NumLoop)*(NumBoard*NumDet*NumSource)*(SubLen+NumBin*DataSize(itry)))
        datatype = DataType{itry};
        datasize = DataSize(itry);
        break;
    end
    if ForceReading == true && isdatatypeargin
        break;
    end
    if ForceReading == false && itry == numel(DataType)
        errordlg({'Can''t handle sizemismatch. Insert more argin' 'Or use (...''forcereading'',''true'') argin'});
        fclose(fid);
        return;
    end
    if ForceReading == true && itry == numel(DataType)
        fh = figure('NumberTitle','off','Name','Choose type','Toolbar','none','menubar','none','HandleVisibility','off','Units','normalized','Position',[0.5 0.5 0.1 0.3]);
        movegui(fh,'center');
        uph = uipanel(fh,'Title','Choose type','units','normalized','position',[0 0 1 1]);
        bg = uibuttongroup(uph,'Visible','on','Position',[0 0 1 1]);
        uicontrol(bg,'style','radiobutton','String','ushort','units','normalized','position',[0 0 1 0.5]);
        uicontrol(bg,'style','radiobutton','String','uint32','units','normalized','position',[0 1/3 1 0.5]);
        uicontrol(bg,'style','radiobutton','String','double','units','normalized','position',[0 2/3 1 0.5]);
        uicontrol(uph,'style','pushbutton','String','Ok','units','normalized','position',[0.5 0 0.5 0.1],'Callback',@AssignDataType);
        waitfor(fh);
        if SkipSub==false
            [NumSource,NumDet,NumBoard]=GetSourceDetBoard(fid,HeadLen,SubLen,NumBin,datatype);
        end
        break;
    end
end

fclose(fid);
NumLoop=CompiledHeader.LoopNum;
CompiledHeader.NumBoard = NumBoard;
CompiledHeader.NumDet = NumDet;
CompiledHeader.NumSource = NumSource;
fid=fopen(FilePath,'rb');
Head=fread(fid,HeadLen,'uint8');

%% Read file
Data=zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1),NumBoard,NumDet,NumSource,NumBin);
Sub=zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1),NumBoard,NumDet,NumSource,SubLen);
if isCompileSubHeader == true
    CompiledSub = CreateDummyStruct(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1),NumBoard,NumDet,NumSource);
end
isbreakcond = false;

try
    for il5= 1:NumLoop(5)
        for il4= 1:NumLoop(4)
            for il3= 1:NumLoop(3)
                for il2=1:NumLoop(2)
                    for il1=1:NumLoop(1)
                        for iB = 1:NumBoard
                            for iD = 1: NumDet
                                for iS = 1:NumSource
                                    if SkipSub == false
                                        BuffSub = fread(fid,SubLen,'uint8');
                                        if ~isempty(BuffSub)&&numel(BuffSub)==SubLen
                                            Sub(il5,il4,il3,il2,il1,iB,iD,iS,:)=BuffSub;
                                        else
                                            warning('backtrace','off')
                                            warnstr = 'Reading interrupted at:\n';
                                            warnstr = strcat(warnstr,strcat('Loop5: ',num2str(il5),'/',num2str(NumLoop(5)),'\n'));
                                            warnstr = strcat(warnstr,strcat('Loop4: ',num2str(il4),'/',num2str(NumLoop(4)),'\n'));
                                            warnstr = strcat(warnstr,strcat('Loop3: ',num2str(il3),'/',num2str(NumLoop(3)),'\n'));
                                            warnstr = strcat(warnstr,strcat('Loop2: ',num2str(il2),'/',num2str(NumLoop(2)),'\n'));
                                            warnstr = strcat(warnstr,strcat('Loop1: ',num2str(il1),'/',num2str(NumLoop(1)),'\n'));
                                            warnstr = strcat(warnstr,strcat('NumBoard: ',num2str(iB),'/',num2str(NumBoard),'\n'));
                                            warnstr = strcat(warnstr,strcat('NumDet: ',num2str(iD),'/',num2str(NumDet),'\n'));
                                            warnstr = strcat(warnstr,strcat('NumSource: ',num2str(iS),'/',num2str(NumSource),'\n'));
                                            warnstr = strcat(warnstr,'Last valid point is the previous iteration');
                                            warnstr = strcat(warnstr,'Output data will have the dimension specified in TRS settings (Header)');
                                            warning(warnstr,'');
                                            warning('backtrace','on')
                                            isbreakcond = true;
                                            break;
                                        end
                                        if isCompileSubHeader  == true
                                            CompiledSub(il5,il4,il3,il2,il1,iB,iD,iS) = FillSub(squeeze(Sub(il5,il4,il3,il2,il1,iB,iD,iS,:)));
                                        end
                                    else
                                        BuffSub = 0;
                                    end
                                    BuffData = fread(fid,NumBin,datatype);
                                    if ~isempty(BuffData)&&numel(BuffData)==NumBin
                                        Sub(il5,il4,il3,il2,il1,iB,iD,iS,:)=BuffSub;
                                    else
                                        warning('backtrace','off')
                                        warnstr = 'Reading interrupted at:\n';
                                        warnstr = strcat(warnstr,strcat('Loop5: ',num2str(il5),'/',num2str(NumLoop(5)),'\n'));
                                        warnstr = strcat(warnstr,strcat('Loop4: ',num2str(il4),'/',num2str(NumLoop(4)),'\n'));
                                        warnstr = strcat(warnstr,strcat('Loop3: ',num2str(il3),'/',num2str(NumLoop(3)),'\n'));
                                        warnstr = strcat(warnstr,strcat('Loop2: ',num2str(il2),'/',num2str(NumLoop(2)),'\n'));
                                        warnstr = strcat(warnstr,strcat('Loop1: ',num2str(il1),'/',num2str(NumLoop(1)),'\n'));
                                        warnstr = strcat(warnstr,strcat('NumBoard: ',num2str(iB),'/',num2str(NumBoard),'\n'));
                                        warnstr = strcat(warnstr,strcat('NumDet: ',num2str(iD),'/',num2str(NumDet),'\n'));
                                        warnstr = strcat(warnstr,strcat('NumSource: ',num2str(iS),'/',num2str(NumSource),'\n'));
                                        warnstr = strcat(warnstr,'Last valid point is the previous iteration');
                                        warnstr = strcat(warnstr,'Output data will have the dimension specified in TRS settings (Header)');
                                        warning(warnstr,'');
                                        warning('backtrace','on')
                                        isbreakcond = true;
                                        break;
                                    end
                                    Data(il5,il4,il3,il2,il1,iB,iD,iS,:)=BuffData;
                                end
                                if isbreakcond, break; end
                            end
                            if isbreakcond, break; end
                        end
                        if isbreakcond, break; end
                    end
                    if isbreakcond, break; end
                end
                if isbreakcond, break; end
            end
            if isbreakcond, break; end
        end
        if isbreakcond, break; end
    end
    fclose(fid);
catch ME
    fclose(fid);
    throw(ME);
end

%% Assign outputs
Data = squeeze(Data);
output{1} = Head;
output{2} = CompiledHeader;
output{3} = squeeze(Sub);
if isCompileSubHeader == false
    output{4} = [];
else
    output{4} = squeeze(CompiledSub);
end
output{5} = Sub;
output{6} = datasize;
output{7} = datatype;
if isCompileSubHeader == false
    output{8} = [];
else
    output{8} = CompiledSub;
end

varargout = output(1:NumArgOut);
end

function TS = CreateDummyStruct(Loop5,Loop4,Loop3,Loop2,Loop1,NumBoard,NumDet,NumSource)
FieldNames = {'Geom','Source','Fiber','Det','Board','Coord','Pad','Xf','Yf','Zf','Rf','Xs','Ys','Zs','Rs','Rho','TimeNom','TimeEff'...
    'n','Loop','Acq','Page','RoiNum','RoiFirst','RoiLast','RoiLambda','RoiPower'};
TS = struct.empty(Loop5,Loop4,Loop3,Loop2,Loop1,NumBoard,NumDet,0);
for i5 = 1:Loop5
    for i4 = 1:Loop4
        for i3 = 1:Loop3
            for i2 = 1:Loop2
                for i1 = 1:Loop1
                    for iB = 1:NumBoard
                        for iD = 1: NumDet
                            for iS = 1:NumSource
                                for ifields = 1:numel(FieldNames)
                                TS(i5,i4,i3,i2,i1,iB,iD,iS).(FieldNames{ifields}) = 0;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

end

function FS = FillSub(Sub,varargin)
% Optimized for speed
FieldNames = {'Geom','Source','Fiber','Det','Board','Coord','Pad','Xf','Yf','Zf','Rf','Xs','Ys','Zs','Rs','Rho','TimeNom','TimeEff'...
    'n','Loop','Acq','Page','RoiNum','RoiFirst','RoiLast','RoiLambda','RoiPower'};
Char = 1; CT = 1;
Long = 4; LT = 2;
Double = 8; DT = 3;
Short = 2; ST = 4;
Unit = 1; %byte
MAX_ROI_SUB = 4;
MAX_LOOP = 3;
FieldSize = [Char,Char,Char,Char,Char,Char,Char,Double,Double,Double,Double,Double,Double,Double,Double,Double,Double,Double,Double,MAX_LOOP*Long,Long,Long,Char,MAX_ROI_SUB*Short,MAX_ROI_SUB*Short,MAX_ROI_SUB*Double,MAX_ROI_SUB*Double]./Unit;
FieldType = [CT,CT,CT,CT,CT,CT,CT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,LT,LT,LT,CT,ST,ST,DT,DT];
iF = 1;
block = 0;
if nargin>1
    FS.Source = Sub(2);
    FS.Det = Sub(4);
    FS.Board = Sub(5);
else
    nSub = numel(Sub);
    while block(end)<nSub
        block = block(end) + (1:FieldSize(iF));
        switch FieldType(iF)
            case CT
                FS.(FieldNames{iF}) = Sub(block);
                if strcmp(FieldNames{iF},'Geom')
                    if FS.(FieldNames{iF}) == 0
                        FS.(FieldNames{iF}) = 'REFL';
                    else
                        FS.(FieldNames{iF}) = 'TRASM';
                    end
                end
                if strcmp(FieldNames{iF},'Coord')
                    if FS.(FieldNames{iF}) == 0
                        FS.(FieldNames{iF}) = 'CART';
                    else
                        FS.(FieldNames{iF}) = 'POLAR';
                    end
                end
            case LT
                FS.(FieldNames{iF}) = cast(typecast(uint8(Sub(block)),'int32'),'double');
            case DT
                FS.(FieldNames{iF}) = typecast(uint8(Sub(block)),'double');
            case ST
                FS.(FieldNames{iF}) = cast(typecast(uint8(Sub(block)),'int16'),'double');
        end
        iF = iF +1;
    end
end
end

function FH = FillHeader(Head)
FieldNames = {'Ver','SubHeader','SubHeaderVer','SizeHeader','SizeSubHeader','SizeData','Kind','Appl'...
    ,'Oma','Date','Time','LoopHome','LoopFirst','LoopLast','LoopDelta','LoopNum','McaChannNum','PageNum'...
    'FrameNum','RamNum','McaTime','McaFactor','MeasNorm','LabelName','LabelContent','Constn','ConstRho'...
    'ConstThick','MammHeader','MammIdxFirst','MammIdxLast','MammIdxTop','MammRateMid','MammRateHigh'};
Char = 1; CT = 1;
Long = 4; LT = 2;
Double = 8; DT = 3;
Short = 2; ST = 4;
Unit = 1; %byte
MAX_LOOP = 3;
LABEL_MAX = 16;
LABEL_NAMELEN = 12;
LABEL_CONTENTLEN = 22;
DX = 2;
FieldSize = [2*Short,Long,Long,Long,Long,Long,Long,Long,Long,11*Char,9*Char,MAX_LOOP*Long,MAX_LOOP*Long,MAX_LOOP*Long,MAX_LOOP*Long,MAX_LOOP*Long,Long,Long,Long,Long,Double,Double,Long,LABEL_NAMELEN*LABEL_MAX*Char,LABEL_CONTENTLEN*LABEL_MAX*Char,Double,Double,Double,Long,DX*Long,DX*Long,DX*Long,DX*Long,DX*Long]./Unit;
FieldType = [ST,LT,LT,LT,LT,LT,LT,LT,LT,CT,CT,LT,LT,LT,LT,LT,LT,LT,LT,LT,DT,DT,LT,CT,CT,DT,DT,DT,LT,LT,LT,LT,LT,LT];
iF = 1;
block = 0;
while block(end)<length(Head)
    block = block(end) + (1:FieldSize(iF));
    RawData = Head(block);
    switch FieldType(iF)
        case CT
            if strcmp(FieldNames{iF},'LabelName')
                Data = reshape(RawData,[12,16])';
                Data = char(Data);
                Data = string(Data);
            end
            if strcmp(FieldNames{iF},'LabelContent')
                Data = reshape(RawData,[22,16])';
                Data = char(Data);
                Data = string(Data);
            end
            if strcmp(FieldNames{iF},'Date')
                Data = datetime(char(RawData'),'InputFormat','MM-dd-yyyy');
            end
            if strcmp(FieldNames{iF},'Time')
                Data = char(RawData');
            end
        case LT
            Data = cast(typecast(uint8(RawData),'int32'),'double');
            if strcmp(FieldNames{iF},'Kind')
                switch Data
                    case 0
                        Data = 'Measure';
                    case 1
                        Data = 'System';
                    case 2
                        Data = 'Simul';
                end
            end
            if strcmp(FieldNames{iF},'Appl')
                switch Data
                    case 0
                        Data = 'Diff';
                    case 1
                        Data = 'Mamm';
                    case 2
                        Data = 'Oxym';
                    case 3
                        Data = 'Fluo';
                    case 4
                        Data = 'Spec';
                end
            end
            if strcmp(FieldNames{iF},'SubHeader')
                switch Data
                    case 0
                        Data = false;
                    case 1
                        Data = true;
                end
            end
            if strcmp(FieldNames{iF},'Oma')
                switch Data
                    case 0
                        Data = false;
                    case 1
                        Data = true;
                end
            end
            if strcmp(FieldNames{iF},'MammHeader')
                switch Data
                    case 0
                        Data = false;
                    case 1
                        Data = true;
                end
            end
        case DT
            Data = typecast(uint8(RawData),'double');
        case ST
            if strcmp(FieldNames{iF},'Ver')
                Data = cast(typecast(uint8(RawData(1:2)),'int16'),'double');
                Data = horzcat(Data,cast(typecast(uint8(RawData(3:4)),'int16'),'double')); %#ok<AGROW>
            else
                Data = cast(typecast(uint8(RawData),'int16'),'double');
            end
    end
    
    FH.(FieldNames{iF}) = Data;
    iF = iF +1;
end

end

function [output]=string2boolean(string)
if strcmp(string,'false')
    output = false;
else
    output = true;
end
end

function AssignDataType(src,~)
ph = src.Parent;
rbh = findobj(ph,'style','radiobutton');
assignin('caller','datatype',rbh(logical([rbh.Value])).String);
fh = ancestor(src,'figure');
close(fh);
end

function [NumSource,NumDet,NumBoard]=GetSourceDetBoard(fid,HeadLen,SubLen,NumBin,DataType)
frewind(fid);
fread(fid,HeadLen,'uint8');
out = false;
BuffParms = []; Parms = [];
while(out==false)
    SubRaw=fread(fid,SubLen,'uint8');
    if isempty(SubRaw), break; end
    try
        CompSub = FillSub(SubRaw,'quick');
    catch ME
        if strcmpi(ME.identifier,'MATLAB:badsubscript')
            out = true;
        else
            rethrow(ME);
        end
    end
    fread(fid,NumBin,DataType);
    ActParms = [CompSub.Source CompSub.Det CompSub.Board];
    if(isequal(BuffParms,ActParms))
        out = true;
    else
        if isempty(BuffParms), BuffParms = ActParms; end
        if isempty(Parms)
            Parms = ActParms;
        else
            Parms(end+1,:) = ActParms; %#ok<AGROW>
        end
    end
end
NumSource = numel(unique(Parms(:,1)));
NumDet = numel(unique(Parms(:,2)));
NumBoard = numel(unique(Parms(:,3)));
end