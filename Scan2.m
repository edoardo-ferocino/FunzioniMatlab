InitScript
[A,~,CH]=DatRead3('PT_DX0001','ForceReading',true);
[NumY,NumX,NumChan,NumBin]=size(A);
SETT = TRSread('C:\Temp\TRS.TRS');
% for iy = NumY:-1:1
%     if(rem(iy,2)==1)
%         A=flip(A,2);
%     end
% end
AcqTime = CH.McaTime;
FromLeft2RightPatienView = true;
AllCounts = sum(A,4);
if FromLeft2RightPatienView == true
    AllCounts = flip(AllCounts);
end
Xv = linspace(-NumX/2,NumX/2,NumX);
figure
for ich = 1 : NumChan
    subplot(4,2,ich);
    imagesc(AllCounts(:,:,ich)./AcqTime);
    %pcolor(-MapCounts(:,:,ich));
    colormap pink, shading interp, axis equal;
    AX = gca;
    AX.Title.String = num2str(ich);
    %   AX.XTickLabel = Xv(AX.XTick);%linspace(StartX,StopX,numel(AX.XTick));
    %AX.YTickLabel = Yv(AX.YTick);%linspace(StartY,StopY,numel(AX.YTick));
end
FH=figure;
AllChan=sum(AllCounts,3);
imagesc(AllChan./AcqTime);
%pcolor(AllChan);
colormap pink, shading interp, axis equal;
SumChan = squeeze(sum(A,3));
if FromLeft2RightPatienView == true
    SumChan = flip(SumChan);
end
PickCurve(FH,SumChan);

figure, imagesc(flip(flip(AllChan,1),2))
colormap pink, shading interp, axis equal;
figure
Waves =[635 680 785 905 930 975 1060];
for iw = 1:numel(Waves)
    WavwMatrix(iw).Data = A(:,:,:,SETT.Roi(iw,2)+1:SETT.Roi(iw,3)+1);
    for ich = 1:NumChan
        WavwMatrix(iw).Chan(ich).Data = WavwMatrix(iw).Data(:,:,ich,:);
    end
    WavwMatrix(iw).SumChanData = squeeze(sum(WavwMatrix(iw).Data,3));
    WavwMatrix(iw).Curves = WavwMatrix(iw).SumChanData;
    WavwMatrix(iw).CountsAllChan = squeeze(sum(WavwMatrix(iw).Curves,3));
    subplot(4,2,iw);
    imagesc(WavwMatrix(iw).CountsAllChan./AcqTime);
    %pcolor(-MapCounts(:,:,ich));
    colormap pink, shading interp, axis equal;
    AX = gca;
    AX.Title.String = num2str(Waves(iw));
    
end
return
for  iy = NumY-1:-1:0
    NewRow = zeros(1,NumX);
    if(rem(iy,2)==1)
        NewRow(1,1:end-2) = AllChan(iy+1,2+1:end);
    else
        NewRow = AllChan(iy+1,1:end);
    end
    NewAllChan(iy+1,:) = NewRow;
    figure(4), imagesc(NewAllChan); colormap pink,shading interp, axis equal;
end
figure(4), imagesc(NewAllChan); colormap pink,shading interp, axis equal;
return
FromLeft2RightPatienView = false;


NumX = 160; StartX = -80; StopX = 80; PaceX = 1; Xv = StartX:PaceX:StopX;
NumY = 20; StartY = 0;StopY = -19; PaceY = -1; Yv = StartY:PaceY:StopY;
Raw = DatRead2('kjdfhdhs0000',1215,NumChan,NumX,NumY,'subheader',false);
if isstruct(Raw), return, end;
h=figure;
semilogy(squeeze(Raw(1,NumX/2,1,:)))
AddPickX(h);
waitfor(h);
PP = Zone(1).PickPoints{1};
warning('off');
Counts = zeros(NumY,NumX);
for iy = 1:NumY
    for ix = 1:NumX
        for ich = 1:NumChan
            Counts(iy,ix) = Counts(iy,ix) + sum(Raw(iy,ix,ich,:));
            Raw2(iy,ix,ich,:) = Raw(iy,ix,ich,:) - mean(Raw(iy,ix,ich,(PP(1):PP(2))));
        end
    end
end
warning('on')
semilogy(squeeze(Raw2(1,NumX/2,1,:)))

MapCounts = sum(Raw2,4);

if FromLeft2RightPatienView == true
    MapCounts = flip(MapCounts,2);
    Counts = flip(Counts,2);
end

% subplot(4,2);
figure
for ich = 1 : NumChan
    subplot(4,2,ich);
    %figure(ich)
    %FFS(num2str(ich))
    %FigFullScreen;
    
    imagesc(-MapCounts(:,:,ich));
    
    %pcolor(-MapCounts(:,:,ich));
    AX = gca;
    AX.Title.String = num2str(ich);
    AX.XTickLabel = Xv(AX.XTick);%linspace(StartX,StopX,numel(AX.XTick));
    AX.YTickLabel = Yv(AX.YTick);%linspace(StartY,StopY,numel(AX.YTick));
    colormap pink, shading interp, axis equal;
end
figure
AllChan = sum(MapCounts./AcqTime,3);
imagesc(flip(AllChan));
%pcolor(AllChan);
colormap pink, shading interp, axis equal;
AX = gca;
AX.XTickLabel = linspace(StartX,StopX,numel(AX.XTick));
AX.YTickLabel = linspace(StartY,StopY,numel(AX.YTick));

fh=figure;
    imagesc(-Counts./AcqTime);
    AX = gca;
    AX.XTickLabel = Xv(AX.XTick);%linspace(StartX,StopX,numel(AX.XTick));
    AX.YTickLabel = Yv(AX.YTick);%linspace(StartY,StopY,numel(AX.YTick));
    colormap pink, shading interp, axis equal;

DefineBorder(fh,PaceY,PaceX);
