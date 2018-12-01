InitScript
[RawData, H, CH] = DatRead3('rebecca0000_summed','datatype','uint32','forcereading',true);
[RawIRF, HIRF, CHIRF] = DatRead3('IRF0000_summed','datatype','uint32');
Sett = TRSread('../Settings/TRS.TRS');
Factor = 0.08245;
Wavelenghts = [630 680 785 905 933 975 1060];
FractFirst = 0.1; FractLast = 0.01;  
GateWidth = 500; GateNum = 10;
nW = numel(Wavelenghts);
Data = RawData;
Data(Data==0)=NaN;
[numY,numX,numbin]=size(RawData);

for inr = 1:nW
    range = Sett.Roi(inr,2)+1:Sett.Roi(inr,3)+1;
    Waves(inr).IRF = RawIRF(range);
    [MaxWaveIRFVal,MaxWaveIRFPos]=max(Waves(inr).IRF);
    Waves(inr).MaxIRFPos =MaxWaveIRFPos;%+Sett.Roi(inr,2);
    Waves(inr).Range = range;
end
for inr = 1:nW
    range = Sett.Roi(inr,2)+1:Sett.Roi(inr,3)+1;
    Waves(inr).Data = Data(:,:,range);
    for it = 1:numel(range)
       SingleBinPage = Waves(inr).Data(:,:,it); SingleBinPage=SingleBinPage(:);
       Waves(inr).Median(it) = median(SingleBinPage,'omitnan');
    end
    Waves(inr).Range = range;
    Waves(inr).FwhmMedian = CalcWidth(Waves(inr).Median,0.5)*Factor;
end
for inr = 1:nW
    [MaxVal,MaxPos] = max(Waves(inr).Data,[],3,'omitnan');
    NormWave = squeeze(Waves(inr).Data./MaxVal);
    Waves(inr).NumGates = GateNum;
    Waves(inr).GateWidth = GateWidth;
    Waves(inr).Wavelenghts = Wavelenghts;
    for iy = 1:numY
        for ix = 1:numX
            if ~isnan(sum(NormWave(iy,ix,:)))
                Waves(inr).FwhmData(iy,ix) = CalcWidth(Waves(inr).Data(iy,ix,:),0.5);
                if Waves(inr).FwhmData(iy,ix)>Waves(inr).FwhmMedian
                    Indxes = find(NormWave(iy,ix,:)>=FractFirst);
                    FirstChannel = Indxes(find(Indxes<=MaxPos(iy,ix),1,'first'));
                    Indxes = find(NormWave(iy,ix,:)<=FractLast);
                    LastChannel = Indxes(find(Indxes>=MaxPos(iy,ix),1,'first'));
                    Waves(inr).CuttedData = zeros(numel(Waves(inr).Range),1);
                    Waves(inr).CuttedData(FirstChannel:LastChannel) = Waves(inr).Data(iy,ix,FirstChannel:LastChannel);
                    [GatedCurve,GateIntervals]=CalcGate(Waves(inr).MaxIRFPos,Waves(inr).CuttedData,GateWidth,GateNum);
                    Waves(inr).GatedMatrix(iy,ix,:) = GatedCurve;
                    %Waves(inr).GateIntervals(iy,ix,:) = GateIntervals;
                else
                    Waves(inr).GatedMatrix(iy,ix,:) = zeros(GateNum,1);
                end
            else
                Waves(inr).GatedMatrix(iy,ix,:) = zeros(GateNum,1);
            end
        end
    end
    
    
end

PlotGate
function [Gate,Limits]=CalcGate(MaxPosIrf,Data,GateWidth,GateNum)
 Factor = 0.08245*1000; Data = squeeze(Data);
 NumTimeBin = numel(Data);
 FirstBin = MaxPosIrf; LastBin = NumTimeBin;
 RealLenghtTimeScale = (LastBin-FirstBin)*Factor;
 if GateNum*GateWidth>RealLenghtTimeScale
     warning('Data not consistent. Data will be cut')
     GateWidth = RealLenghtTimeScale/GateNum;
 end
 RealData = Data(FirstBin:LastBin);
 TimeBinGate = floor(GateWidth/Factor);
 for ig = 1:GateNum
    Interval = (1:TimeBinGate)+(ig-1)*TimeBinGate;
    Limits(ig,1) = Interval(1); Limits(ig,2) = Interval(end);
    Gate(ig)=sum(RealData((1:TimeBinGate)+(ig-1)*TimeBinGate),'omitnan');
 end
 Limits = Limits+FirstBin-1;
end