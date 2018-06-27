function CompiledHeader = CompileHeader(CH)
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

CompiledHeader = [];
nfields = numel(FieldNames);
for iF = 1:nfields
    RawData=CH.(FieldNames{iF});
    %isCHAR = false;
    switch FieldType(iF)
        case CT
            if strcmp(FieldNames{iF},'LabelName')
                Data = char(RawData); Data = Data'; Data = Data(:);
                Data=typecast(uint8(Data),'uint8');
                %isCHAR = true;
            end
            if strcmp(FieldNames{iF},'LabelContent')
                Data = char(RawData); Data = Data'; Data = Data(:);
                Data=typecast(uint8(Data),'uint8');
                %isCHAR = true;
            end
            if strcmp(FieldNames{iF},'Date')
                Data = typecast(uint8(RawData),'uint8');
                %isCHAR = true;
            end
            if strcmp(FieldNames{iF},'Time')
                Data = typecast(uint8(RawData),'uint8');
                %isCHAR = true;
            end
            %isCHAR = false;
        case LT
            if strcmp(FieldNames{iF},'Kind')
                switch RawData
                    case 'Measure'
                        RawData = 0;
                    case 'System'
                        RawData = 1;
                    case 'Simul'
                        RawData = 2;
                end
            end
            if strcmp(FieldNames{iF},'Appl')
                switch RawData
                    case 'Diff'
                        RawData = 0;
                    case 'Mamm'
                        RawData = 1;
                    case 'Oxym'
                        RawData = 2;
                    case 'Fluo'
                        RawData = 3;
                    case 'Spec'
                        RawData = 4;
                end
            end
            if strcmp(FieldNames{iF},'Oma') || strcmp(FieldNames{iF},'MammHeader') || strcmp(FieldNames{iF},'SubHeader')
                switch RawData
                    case false
                        RawData = 0;
                    case true
                        RawData = 1;
                end
            end
            if contains(FieldNames{iF},'Loop')
                RawData = flip(RawData(1:3));
            end
            Data = typecast(int32(RawData),'uint8');
        case DT
            Data = typecast(double(RawData),'uint8');
        case ST
            if strcmp(FieldNames{iF},'Ver')
                Data = typecast(int16(RawData(1)),'uint8');
                Data = [Data typecast(int16(RawData(2)),'uint8')];
            else
                Data = typecast(int16(RawData),'uint8');
            end
    end
    [nr,nc]=size(Data); if nr>nc, Data = Data'; end
    CompiledHeader = [CompiledHeader; Data'];
end
CompiledHeader=cast(CompiledHeader,'double');
end