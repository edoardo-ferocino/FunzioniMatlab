function CompiledSubHeader = CompileSubHeader(SH,SubHead)
FieldNames = {'Geom','Source','Fiber','Det','Board','Coord','Pad','Xf','Yf','Zf','Rf','Xs','Ys','Zs','Rs','Rho','TimeNom','TimeEff'...
    'n','Loop','Acq','Page','RoiNum','RoiFirst','RoiLast','RoiLambda','RoiPower'};
C = 1; CT = 1;
L = 4; LT = 2;
D = 8; DT = 3;
S = 2; ST = 4;
Unit = 1; %byte
FieldSize = [C,C,C,C,C,C,C,D,D,D,D,D,D,D,D,D,D,D,D,3*L,L,L,C,4*S,4*S,4*D,4*D]./Unit;
FieldType = [CT,CT,CT,CT,CT,CT,CT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,LT,LT,LT,CT,ST,ST,DT,DT];

nfields = numel(FieldNames);
[Nl1, Nl2, Nl3, Nl4, Nl5, Nb, Nd] = size(SH);
for il1 = 1:Nl1
    for il2 = 1:Nl2
        for il3 = 1:Nl3
            for il4 = 1:Nl4
                for il5 = 1:Nl5
                    for ib = 1:Nb
                        for id = 1:Nd
                            CompiledSubHeaderBuffer = [];
                            RD = SH(il1,il2,il3,il4,il5,ib,id);
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
                                
                                [nr,nc]=size(Data); if nr>nc, Data = Data'; end
                                CompiledSubHeaderBuffer = [CompiledSubHeaderBuffer; Data'];
                                if ~all(CompiledSubHeaderBuffer == squeeze(SubHead(il1,il2,il3,il4,il5,ib,id,1:sum(FieldSize(1:iF)))))
                                    kk = 2;
                                end
                            end
                            CompiledSubHeaderBuffer=cast(CompiledSubHeaderBuffer,'double');
                            CompiledSubHeader(il1,il2,il3,il4,il5,ib,id,:) = CompiledSubHeaderBuffer;
                        end
                    end
                end
            end
        end
    end
end
end