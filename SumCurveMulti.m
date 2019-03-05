InitScript

%%
IRF_FN = 'IRF0000';
Data_FN = 'rebecca0001';
is_scan = true; is_sum_rep_irf = true; is_sum_rep_data = true;
SETT = TRSread('../Settings/TRS');
[IRF,H,CH,SH,CSH,UnSquezSubs]=DatRead3(IRF_FN,'compilesub',true);
ndIRF = ndims(IRF);
if ndIRF == 2
    [NumChan,NumBin]=size(IRF);
    NumRep = 1;
elseif ndIRF == 3
    [NumRep,NumChan,NumBin] = size(IRF);
    if is_sum_rep_irf
        IRF = (sum(IRF,1));
    end
else
    return
end

POS_IRF_raw = zeros(NumChan,1);
for ich = 1:NumChan
    StartBin = SETT.Roi(1,2)+1;
    StopBin = SETT.Roi(1,3)+1;
    B=squeeze(IRF(1,ich,StartBin:StopBin));
    [width,~,maxpos] = CalcWidth(B,0.5);
    pos = maxpos;
    pos = pos + StartBin-1;
    %pos = pos + round((iw-1)*Period/Factor);
    POS_IRF_raw(ich) = pos;
end

mean_pos = mean(POS_IRF_raw);
shift = mean_pos - POS_IRF_raw;
IRF_Shifted = zeros(size(IRF));
for ich = 1:NumChan
    IRF_Shifted(:,ich,:)=circshift(IRF(:,ich,:),round(shift(ich)),3);
end

IRF_Shifted_summed = (sum(IRF_Shifted,2));
save('shift.mat','shift');

fid_out = fopen([IRF_FN '_summed.DAT'], 'wb');
if is_sum_rep_irf
    CH.LoopNum(1) = 1; CH.LoopLast(1) = 1;
    H=CompileHeader(CH);
    NumRep = 1;
end
fwrite(fid_out, H, 'uint8');
for irep = 1:NumRep
    fwrite(fid_out, SH(irep,1,:), 'uint8');
    curve=squeeze(IRF_Shifted_summed(irep,1,:));
    fwrite(fid_out, curve, 'uint32');
end
warning('\nConversion to long (uint32) required for file: %s', [IRF_FN '_summed']);
fclose(fid_out);


%%
[Data,H,CH,SH,CSH] = DatRead3(Data_FN,'forcereading',true);
if is_scan == 1
    [NumY,NumX,NumChan,NumBin]=size(Data);
    Data_Shifted = zeros(NumY,NumX,NumChan,NumBin);
    for iy =1: NumY
        for ix = 1: NumX
            for ich=1:NumChan
                if sum(Data(iy,ix,ich,:))~=0
                    Data_Shifted(iy,ix,ich,:)=circshift(squeeze(Data(iy,ix,ich,:)),round(shift(ich)));
                end
            end
        end
    end
else
    ndData=ndims(Data);
    if ndData == 3
        NumRep = size(Data,1);
        if is_sum_rep_data
            Data = (sum(Data,1));
            NumRep = 1;
        end
    end
    Data_Shifted = zeros(size(Data));
    for ich=1:NumChan
        Data_Shifted(:,ich,:)=circshift(Data(:,ich,:),round(shift(ich)),3);
    end
end

if is_scan
    Data_Shifted_Summed = squeeze(sum(Data_Shifted,3));
else
    Data_Shifted_Summed = (sum(Data_Shifted,2));
end

fid_out = fopen([Data_FN '_summed.DAT'], 'wb');

if is_scan == 1
    fwrite(fid_out, H, 'uint8');
    for iy = 1:NumY
        for ix = 1:NumX
            fwrite(fid_out, SH(iy,ix,1,:), 'uint8');
            curve=Data_Shifted_Summed(iy,ix,:);
            fwrite(fid_out, curve, 'uint32');
        end
        warning('\nConversion to long (uint32) required for file: %s', [Data_FN '_summed']);
    end
else
    if is_sum_rep_data
        CH.LoopNum(1) = 1; CH.LoopLast(1) = 1;
        H=CompileHeader(CH);
        NumRep = 1;
    end
    fwrite(fid_out, H, 'uint8');
    for irep = 1:NumRep
        fwrite(fid_out, SH(irep,1,:), 'uint8');
        curve=Data_Shifted_Summed(irep,1,:);
        fwrite(fid_out, curve, 'uint32');
    end
end

fclose(fid_out);

