function [ Data, varargout ] = DatRead3(FileName,varargin)
%DatRead3('FileName')
%Can be as input a selection of the following parameters
%DatRead3(...,'loop4',numloop4,'loop5',numloop5,'datatype','uint32','compilesub',true/false,'forcereading',true/false)
%[Data,Header,EasyReadableHead,SubHeaders,EasyReadableSubHead,UnSqueezedHeader]=DatRead3(...)

NumArgOut = nargout-1;
NumArgin = nargin-1;
sumloop = 0; ForceReading = false;

isCompileSubHeader = false;
HeadLen=764;
FilePath = [FileName '.DAT'];
if isempty(fileparts(FileName))
    FilePath = fullfile(pwd,[FileName,'.DAT']);
end
fid=fopen(FilePath,'rb');
if fid<0, errordlg('File not found'); Data = []; return; end %#ok<UNRCH>
Head=fread(fid,HeadLen,'uint8');
datatype = 'ushort';
CompiledHeader = FillHeader(Head);
SubLen=CompiledHeader.SizeSubHeader;
NumBin = CompiledHeader.McaChannNum;
CompiledHeader.LoopNum(4) = 1; CompiledHeader.LoopNum(5) = 1;
CompiledHeader.LoopFirst(4) = 0; CompiledHeader.LoopFirst(5) = 0;
CompiledHeader.LoopLast(4) = 0; CompiledHeader.LoopLast(5) = 0;
CompiledHeader.LoopDelta(4) = 1; CompiledHeader.LoopDelta(5) = 1;
CompiledHeader.LoopHome(4) = 0; CompiledHeader.LoopHome(5) = 0;
CompiledHeader.LoopNum(1:3)=flip(CompiledHeader.LoopNum(1:3));
CompiledHeader.LoopFirst(1:3)=flip(CompiledHeader.LoopFirst(1:3));
CompiledHeader.LoopLast(1:3)=flip(CompiledHeader.LoopLast(1:3));
CompiledHeader.LoopDelta(1:3)=flip(CompiledHeader.LoopDelta(1:3));
CompiledHeader.LoopHome(1:3)=flip(CompiledHeader.LoopHome(1:3));

for iN = 1:NumArgin
    if strcmpi(varargin{iN},'loop5')
        CompiledHeader.LoopNum(5) = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'loop4')
        CompiledHeader.LoopNum(4) = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'loop3')
        CompiledHeader.LoopNum(3) = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'loop2')
        CompiledHeader.LoopNum(2) = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'loop1')
        CompiledHeader.LoopNum(1) = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'datatype')
        datatype = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'compilesub')
        isCompileSubHeader = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'sumoverloop')
        sumloop = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'forcereading')
        ForceReading = varargin{iN+1};
    end
end
if strcmpi(datatype,'uint32')
    datatry = {'uint32'};
else
    datatry = {'ushort','uint32'};
end
for itry = 1:numel(datatry)
    SubRaw=fread(fid,SubLen,'uint8'); CompSub = FillSub(SubRaw); fread(fid,NumBin,datatry{itry});
    BuffBoard = CompSub.Board; BuffDet = CompSub.Det; nDet = 1; nBoard = 1;
    out = false;
    while(out==false)
        SubRaw=fread(fid,SubLen,'uint8');
        if isempty(SubRaw), break; end
        CompSub = FillSub(SubRaw); fread(fid,NumBin,datatype);
        if(BuffBoard == CompSub.Board && BuffDet == CompSub.Det)
            out = true;
        else
            if BuffDet ~= CompSub.Det, nDet = nDet+1; end
            if BuffBoard ~= CompSub.Board, nBoard = nBoard+1; end
        end
    end
    
    NumLoop=CompiledHeader.LoopNum;
    
    info=dir(FilePath);
    if info.bytes ~= (HeadLen + prod(NumLoop)*(nBoard*nDet)*(SubLen+NumBin*2))
        nBoardBuff = nBoard; nDetBuff = nDet;
        if info.bytes == (HeadLen + prod(NumLoop)*(nBoard*nDet)*(SubLen+NumBin*4))
            if ~strcmpi(datatype,'uint32')
                datatype = 'uint32';
                warning('Datasize was ''uint32''. Data will be read. To be sure, use the argument ''datatype'' and set it to ''uint32''');
                %fclose(fid);
                %return
            end
        else
            if (ForceReading==false)&&itry==numel(datatry)
                nBoard = nBoardBuff; nDet = nDetBuff; datatype = 'ushort';
                RemBytes = info.bytes-(HeadLen + prod(NumLoop)*(nBoard*nDet)*(SubLen+NumBin*2));
                RemLoop = RemBytes/(prod(NumLoop)*(nBoard*nDet)*(SubLen+NumBin*2));
                CompiledHeader.LoopNum(4) = RemLoop+1; NumLoop(4) = RemLoop+1;
                NumLoop(5) = 1;
                warndlg({['For the 4 loop the value ',num2str(NumLoop(4)),' will be used. Loop5 is undefined'],....
                    'If you want to insert the correct values for loop4 and loop5 launch again the function with that arguments'},'Dimension mismatch','modal')
            end
            if (ForceReading==true)&&itry==numel(datatry)
                prompt = {'Enter datatype'};
                title = 'Datatype';
                dims = [1 35];
                definput = {'ushort/uint32'};
                answer = inputdlg(prompt,title,dims,definput);
                datatype = answer{1}; 
                fclose(fid); fid=fopen(FilePath,'rb');
                Head=fread(fid,HeadLen,'uint8');
                SubRaw=fread(fid,SubLen,'uint8'); CompSub = FillSub(SubRaw); fread(fid,NumBin,datatype);
                BuffBoard = CompSub.Board; BuffDet = CompSub.Det; nDet = 1; nBoard = 1;
                out = false;
                while(out==false)
                    SubRaw=fread(fid,SubLen,'uint8');
                    if isempty(SubRaw), break; end
                    CompSub = FillSub(SubRaw); fread(fid,NumBin,datatype);
                    if(BuffBoard == CompSub.Board && BuffDet == CompSub.Det)
                        out = true;
                    else
                        if BuffDet ~= CompSub.Det, nDet = nDet+1; end
                        if BuffBoard ~= CompSub.Board, nBoard = nBoard+1; end
                    end
                end
                NumLoop=CompiledHeader.LoopNum;
            end
        end
        frewind(fid);
        Head=fread(fid,HeadLen,'uint8');
    else
        break
    end
end
fclose(fid);


fid=fopen(FilePath,'rb');
Head=fread(fid,HeadLen,'uint8');

A=zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1),nBoard,nDet,NumBin);
Sub=zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1),nBoard,nDet,SubLen);
%CompiledSub = zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1));
try
    for il5= 1:NumLoop(5)
        for il4= 1:NumLoop(4)
            for il3= 1:NumLoop(3)
                for il2=1:NumLoop(2)
                    for il1=1:NumLoop(1)
                        for iB = 1:nBoard
                            for iD = 1: nDet
                                TrashSub = fread(fid,SubLen,'uint8');
                                if ~isempty(TrashSub)
                                    Sub(il5,il4,il3,il2,il1,iB,iD,:)=TrashSub;
                                else
                                    warning('backtrace','off')
                                    warning('Reading interrupted at:');
                                    warning(strcat('Loop5: ',num2str(il5),'/',num2str(NumLoop(5))));
                                    warning(strcat('Loop4: ',num2str(il4),'/',num2str(NumLoop(4))));
                                    warning(strcat('Loop3: ',num2str(il3),'/',num2str(NumLoop(3))));
                                    warning(strcat('Loop2: ',num2str(il2),'/',num2str(NumLoop(2))));
                                    warning(strcat('Loop1: ',num2str(il1),'/',num2str(NumLoop(1))));
                                    warning(strcat('NumBoard: ',num2str(iB),'/',num2str(nBoard)));
                                    warning(strcat('NumDet: ',num2str(iD),'/',num2str(nDet)));
                                    warning('Output data will have the dimension specified in TRS settings (Header)');
                                    warning('backtrace','on')
                                    break;
                                end
                                if (isCompileSubHeader)
                                    CompiledSub(il5,il4,il3,il2,il1,iB,iD) = FillSub(squeeze(Sub(il5,il4,il3,il2,il1,iB,iD,:)));
                                end
                                A(il5,il4,il3,il2,il1,iB,iD,:)=fread(fid,NumBin,datatype);
                            end
                            if isempty(TrashSub), break; end
                        end
                        if isempty(TrashSub), break; end
                    end
                    if isempty(TrashSub), break; end
                end
                if isempty(TrashSub), break; end
            end
            if isempty(TrashSub), break; end
        end
        if isempty(TrashSub), break; end
    end
    
    CompiledHeader.NumBoard = nBoard;
    CompiledHeader.NumDet = nDet;
    
    if sumloop
        A=sum(squeeze(A),sumloop);
    end
    
    fclose(fid);
catch ME
    fclose(fid);
    %ME = MException('FileReading:GeneralError',{'Error reading the file:',ME.message, ' Encountered error at Loop1: ', num2str(il1), ' Loop2: ', num2str(il2), ' Loop3: ', num2str(il3), 'Loop4: ', num2str(il4),' Loop5: ', num2str(il5)});
    throw(ME);
    return
end

Data = squeeze(A);
switch NumArgOut
    case 1
        varargout{1}= Head;
    case 2
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
    case 3
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
        varargout{3} = squeeze(Sub);
    case 4
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
        varargout{3} = squeeze(Sub);
        if isCompileSubHeader == false
            varargout{4} = [];
        else
            varargout{4} = squeeze(CompiledSub);
        end
    case 5
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
        varargout{3} = squeeze(Sub);
        if isCompileSubHeader == false
            varargout{4} = [];
        else
            varargout{4} = squeeze(CompiledSub);
        end
        varargout{5} = Sub;
    case 6
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
        varargout{3} = squeeze(Sub);
        if isCompileSubHeader == false
            varargout{4} = [];
        else
            varargout{4} = squeeze(CompiledSub);
        end
        varargout{5} = Sub;
        varargout{6} = datatype;
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
            Data = cast(typecast(uint8(RawData),'int32'),'double');
        case DT
            Data = typecast(uint8(RawData),'double');
        case ST
            Data = cast(typecast(uint8(RawData),'int16'),'double');
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
                Data = [Data cast(typecast(uint8(RawData(3:4)),'int16'),'double')];
            else
                Data = cast(typecast(uint8(RawData),'int16'),'double');
            end
    end
    
    FH.(FieldNames{iF}) = Data;
    iF = iF +1;
end

end
