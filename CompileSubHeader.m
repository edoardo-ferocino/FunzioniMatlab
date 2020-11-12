function CompiledSubHeader = CompileSubHeader(SH)
FieldNames = {'Geom','Source','Fiber','Det','Board','Coord','Pad','Xf','Yf','Zf','Rf','Xs','Ys','Zs','Rs','Rho','TimeNom','TimeEff'...
    'n','Loop','Acq','Page','RoiNum','RoiFirst','RoiLast','RoiLambda','RoiPower'};
CT = 1;
LT = 2;
DT = 3;
ST = 4;
FieldType = [CT,CT,CT,CT,CT,CT,CT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,LT,LT,LT,CT,ST,ST,DT,DT];

nfields = numel(FieldNames);

[Nl1, Nl2, Nl3, Nl4, Nl5, Nb, Nd, Ns] = size(SH);
CompiledSubHeader = zeros([Nl1, Nl2, Nl3, Nl4, Nl5, Nb, Nd, Ns, 204]);

for il1 = 1:Nl1
    for il2 = 1:Nl2
        for il3 = 1:Nl3
            for il4 = 1:Nl4
                for il5 = 1:Nl5
                    for ib = 1:Nb
                        for id = 1:Nd
                            for is = 1:Ns
                                nstart = 0;
                                CompiledSubHeaderBuffer = zeros(204,1);
                                RD = SH(il1,il2,il3,il4,il5,ib,id,is);
                                for iF = 1:nfields
                                    RawData=RD.(FieldNames{iF});
                                    switch FieldType(iF)
                                        case CT
                                            if strcmp(FieldNames{iF},'Geom')
                                                if strcmpi(RawData,'REFL')
                                                    RawData = 0;
                                                elseif strcmpi(RawData,'TRASM')
                                                    RawData = 1;
                                                end
                                            end
                                            if strcmp(FieldNames{iF},'Coord')
                                                if strcmpi(RawData,'CART')
                                                    RawData = 0;
                                                elseif strcmpi(RawData,'POLAR')
                                                    RawData = 1;
                                                end
                                            end
                                            Data = typecast(uint8(RawData),'uint8');
                                        case LT
                                            Data = typecast(int32(RawData),'uint8');
                                        case DT
                                            Data = typecast(double(RawData),'uint8');
                                        case ST
                                            Data = typecast(int16(RawData),'uint8');
                                    end
                                    
                                    if isrow(Data)
                                        Data = Data';
                                    end
                                    CompiledSubHeaderBuffer(nstart+(1:numel(Data)),1) = Data;
                                    nstart = nstart + numel(Data);
                                end
                                CompiledSubHeaderBuffer=cast(CompiledSubHeaderBuffer,'double');
                                CompiledSubHeader(il1,il2,il3,il4,il5,ib,id,is,:) = CompiledSubHeaderBuffer;
                            end
                        end
                    end
                end
            end
        end
    end
end
CompiledSubHeader = squeeze(CompiledSubHeader);
end