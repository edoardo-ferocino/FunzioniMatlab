function [ Data, varargout ] = DatRead2(FileName,NumBin,NumLoop1,varargin)
% e.g. 4096 bin, 8 channel measure, 5 repetisions
% [Data]=DatRead2('Exemple',4096,8,5). If you want to sum the repetitions
% do: [Data]=DatRead2('Exemple',4096,8,5,'sum',2)
%[Data]=DatRead2(FileName,NumBin,NumLoop1,...,NumLoop5)
%[Data]=DatRead2(...,'Sum',LoopIndex) LoopIndex is 1 for Loop1,2 for Loop2
%etc
%[Data]=DatRead2(...,'All') Default is All
%[Data,Header,SubHeaders,EasyReadableHead,EasyReadableSubHead]=DatRead2(...)

NumArgIn = length(varargin);
NumArgOut = nargout-1;
if iscellstr(FileName)
    NumFile = length(FileName);
    BufferData = cell(NumFile,1);
    BufferHead = cell(NumFile,1);
    BufferSub = cell(NumFile,1);
    BufferCompiledSub = cell(NumFile,1);
    BufferCompiledHeader = cell(NumFile,1);
else
    NumFile = 1;
end
for ifile=1:NumFile
HeadLen=764;
SubLen=204;
if NumFile>1
    FilePath=[FileName{ifile} '.DAT'];
else
    FilePath=[FileName '.DAT'];
end
fid=fopen(FilePath,'rb');
if fid<0, error('File not found'); return; end %#ok<UNRCH>
Head=fread(fid,HeadLen,'uint8');
CompiledHeader = FillHeader(Head);

NumLoop = ones(5,1);
NumLoop(1)=NumLoop1;
isCompileSubHeader = true;
datatype = 'ushort';
SumArg.Yes=0;
SumArg.LoopIndex=0;
il=2;
for in=1:NumArgIn
    if(strcmpi(varargin{in},'uint32'))
        datatype = varargin{in}; 
    end
    if(strcmpi(varargin{in},'subheader'))
        isCompileSubHeader = varargin{in+1}; 
    end
    if(strcmpi(varargin{in},'sum'))
        SumArg.Yes = in+1;
        SumArg.LoopIndex =6-varargin{in+1};
    end
    if(isnumeric(varargin{in})&&SumArg.Yes~=in)
        NumLoop(il)=varargin{in};
        il=il+1;
    end
end

A=zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1),NumBin);
Sub=zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1),SubLen);
%CompiledSub = zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1));
% info=dir(FilePath);
% if(info.bytes*8==((NumRep*NumChan)*(16*NumBin+8*SubLen)+HeadLen*8))
%     datatype = 'ushort';
% else
%     datatype = 'uint32';
% end
try
for il5= 1:NumLoop(5)
    for il4= 1:NumLoop(4)
        for il3= 1:NumLoop(3)
            for il2=1:NumLoop(2)
                for il1=1:NumLoop(1)
                    if CompiledHeader.SubHeader
                        Sub(il5,il4,il3,il2,il1,:)=fread(fid,SubLen,'uint8');
                        if (isCompileSubHeader), CompiledSub(il5,il4,il3,il2,il1) = FillSub(squeeze(Sub(il5,il4,il3,il2,il1,:))); end
                    end
                    A(il5,il4,il3,il2,il1,:)=fread(fid,NumBin,datatype);
                end
            end
        end
    end
end
fclose(fid);
catch ME
    fclose(fid);
    errordlg({ME.message, ' Encountered error at Loop1: ', num2str(il1), ' Loop2: ', num2str(il2), ' Loop3: ', num2str(il3), 'Loop4: ', num2str(il4),' Loop5: ', num2str(il5)},'Error','modal');
    Data.Message = {ME.message, ' Encountered error at Loop1: ', num2str(il1), ' Loop2: ', num2str(il2), ' Loop3: ', num2str(il3), 'Loop4: ', num2str(il4),' Loop5: ', num2str(il5)};
    return
end
if(SumArg.Yes)
    Data=squeeze(sum(A,SumArg.LoopIndex)); 
else
    Data = squeeze(A); 
end
if ifile>1, BufferData{ifile} = Data; BufferHead{ifile} = Head; if (CompiledHeader.SubHeader), BufferSub{ifile} = Sub; if (isCompileSubHeader), BufferCompiledSub{ifile} = CompiledSub; end; end; BufferCompiledHeader{ifile} = CompiledHeader; end
end
if ifile>1
    Data = BufferData;
    Head = BufferHead;
    if (CompiledHeader.SubHeader)
        Sub = BufferSub;
        if (isCompileSubHeader), CompiledSub = BufferCompiledSub; end
    end
end
switch NumArgOut
    case 1
        varargout{1}= Head;
    case 2
        varargout{1} = Head;
        if (CompiledHeader.SubHeader), varargout{2} = squeeze(Sub); end
    case 3
        varargout{1} = Head;
        if (CompiledHeader.SubHeader), varargout{2} = squeeze(Sub); end
        varargout{3} = CompiledHeader;
    case 4
        varargout{1} = Head;
        if (CompiledHeader.SubHeader), varargout{2} = squeeze(Sub); end
        varargout{3} = CompiledHeader;
        if (CompiledHeader.SubHeader), if(isCompileSubHeader), varargout{4} = squeeze(CompiledSub); end; end;
end
end

function FS = FillSub(Sub)
FieldNames = {'Geom','Source','Fiber','Det','Board','Coord','Pad','Xf','Yf','Zf','Rf','Xs','Ys','Zs','Rs','Rho','TimeNom','TimeEff'...
    'n','Loop','Acq','Page','RoiNum','RoiFirst','RoiLast','RoiLambda','RoiPower'};
C = 1; CT = 1;
L = 4; LT = 2;
D = 8; DT = 3;
S = 2; ST = 4;
Unit = 1; %byte
FieldSize = [C,C,C,C,C,C,C,D,D,D,D,D,D,D,D,D,D,D,D,3*L,L,L,C,4*S,4*S,4*D,4*D]./Unit;
FieldType = [CT,CT,CT,CT,CT,CT,CT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,LT,LT,LT,CT,ST,ST,DT,DT];
iS = 1;
iF = 1;

while iS<=length(Sub)
iStart = iS;
iStop = iS + FieldSize(iF);
iStop = iStop-1;
iS = iStop+1;
RawData = Sub(iStart:iStop);

switch FieldType(iF)
    case CT
        Data = RawData;
        if strcmp(FieldNames{iF},'Geom')
            if Data == 0
                Data = 'REFL';
            else
                Data = 'TRASM';
            end
        end
        if strcmp(FieldNames{iF},'Coord')
            if Data == 0
                Data = 'CART';
            else
                Data = 'POLAR';
            end
        end
    case LT
        Data = cast(typecast(uint8(RawData),'uint32'),'double');
    case DT
        Data = typecast(uint8(RawData),'double');
    case ST
        Data = cast(typecast(uint8(RawData),'uint16'),'double');
end

FS.(FieldNames{iF}) = Data;
iF = iF +1;
end

end


function FH = FillHeader(Head)
FieldNames = {'Ver','SubHeader','SubHeaderVer','SizeHeader','SizeSubHeader','SizeData','Kind','Appl'...
    ,'Oma','Date','Time','LoopHome','LoopFirst','LoopLast','LoopDelta','LoopNum','McaChannNum','PageNum'...
    'FrameNum','RamNum','McaTime','McaFactor','MeasNorm','LabelName','LabelContent','Constn','ConstRho'...
    'ConstThick','MammHeader','MammIdxFirst','MammIdxLast','MammIdxTop','MammRateMid','MammRateHigh'};
C = 1; CT = 1;
L = 4; LT = 2;
D = 8; DT = 3;
S = 2; ST = 4;
Unit = 1; %byte
FieldSize = [2*S,L,L,L,L,L,L,L,L,11*C,9*C,3*L,3*L,3*L,3*L,3*L,L,L,L,L,D,D,L,192*C,352*C,D,D,D,L,2*L,2*L,2*L,2*L,2*L]./Unit;
FieldType = [ST,LT,LT,LT,LT,LT,LT,LT,LT,CT,CT,LT,LT,LT,LT,LT,LT,LT,LT,LT,DT,DT,LT,CT,CT,DT,DT,DT,LT,LT,LT,LT,LT,LT];
iS = 1;
iF = 1;

while iS<=length(Head)
iStart = iS;
iStop = iS + FieldSize(iF);
iStop = iStop-1;
iS = iStop+1;
RawData = Head(iStart:iStop);
%isCHAR = false;
switch FieldType(iF)
    case CT
        if strcmp(FieldNames{iF},'LabelName')
            Data = reshape(RawData,[12,16])';
            Data = char(Data);
            Data = string(Data);
            %isCHAR = true;
        end
        if strcmp(FieldNames{iF},'LabelContent')
            Data = reshape(RawData,[22,16])';
            Data = char(Data);
            Data = string(Data);
            %isCHAR = true;
        end
        if strcmp(FieldNames{iF},'Date')
            Data = char(RawData');
            %isCHAR = true;
        end
        if strcmp(FieldNames{iF},'Time')
            Data = char(RawData');
            %isCHAR = true;
        end
            %isCHAR = false;
    case LT
        Data = cast(typecast(uint8(RawData),'uint32'),'double');
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
            Data = [Data cast(typecast(uint8(RawData(3:4)),'uint16'),'double')];
        else
            Data = cast(typecast(uint8(RawData),'uint16'),'double');
        end
end

FH.(FieldNames{iF}) = Data;
iF = iF +1;
end

end
