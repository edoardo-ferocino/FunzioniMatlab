%% Init
InitScript
% FileNames = {'Irft0043','Irfm0043','Irft0044','Irfm0044','Irft0045','Irfm0045'};
% SeqFileNames = {'Irf0043_seqt','Irf0043_seqm','Irf0044_seqt','Irf0044_seqm','Irf0045_seqt','Irf0045_seqm'};
% AreaTableFileNames = {'area_map_pinhole.txt','area_map_pinhole.txt','area_map_pinhole.txt','area_map_pinhole.txt','area_map_pinhole.txt','area_map_pinhole.txt'};
FileNames = {'Irft0049','Pham0002','Pham0003','Pham0004'};
SeqFileNames = {'Irf0049_seqt','Pha0002_seqm','Pha0003_seqm','Pha0004_seqm'};
AreaTableFileNames = {'area_map_pinhole.txt','area_map_pinhole.txt','area_map_pinhole.txt','area_map_pinhole.txt','area_map_pinhole.txt','area_map_pinhole.txt'};
DNLCoeffFileName = 'coeff_DNL.mat';
dirdata = strsplit(mfilename('fullpath'),'\'); dirdata = dirdata{end-1};
dirdata = strcat('..\..\Data\',dirdata,'\');
dirsett = '..\..\Settings\Recon\';
LITE = 1; MEDIUM = 2; HARD = 3;

%% Options
IsWriteCurveToFile = false;
LoadData = true;
CompensateCountingLoss = true;
DNLCorrection = true;
RemoveBkg = true;
AutoMask = true;
isPlot = false;
LevelIsPlot = MEDIUM;

%% Variables
DelayRange = Inf;
MaskWidthPs = 200;%maggiore dello step tra gate per avere overlap
MaskStart_ch = [45, 45, 45, 50]; % uno per ogni lambda. canale da cui iniziare a tagliare
BkgWeigth = 1;
PercMaxAutoMask = 0.9;
AP_prob = 0;
% Only the last 1/3 of the repetitions are taken as valid
% Factor is set to 80 ps/ch!!!!
%% Constants
FirstValidChan = 1;
FinalValidChan = 120;
FinalDummyLineID = 42;
DummyLineID = 21;
OPTODE = 2;
AreaCh = 125;
IntensityCh = [128 127];
TempChan = 121;
NumChanToPad = 256;
DelayCourse = [0:7; 0:500:500*7]';
DelayFine = [0 0;68 97.5;100 200;128 308;156 397];
col = ['g','b','r','g','k','b','r','g'];
ValidChanRange = FirstValidChan:FinalValidChan;
NumValidChan = numel(ValidChanRange);
warning('off','MATLAB:Axes:NegativeDataInLogAxis');

for ifn = 1:numel(FileNames)
    %% Paths
    FileName = FileNames{ifn};
    SeqFileName = SeqFileNames{ifn};
    AreaTableFileName = AreaTableFileNames{ifn};
    FilePath = strcat(dirdata,FileName,'.dat');
    MatFilePath = strcat(dirdata,FileName,'.mat');
    SeqFilePath = strcat(dirsett,SeqFileName,'.ini');
    SpadMeasFilePath = strcat(dirdata,AreaTableFileName);
    
    %% Read file
    SPAD_meas = readmatrix(SpadMeasFilePath,'Delimiter','\t','FileType','text');
    SPAD_meas = SPAD_meas(:,5);
    Seq = readmatrix(SeqFilePath,'Delimiter','\t','FileType','text');
    load(strcat(dirdata,DNLCoeffFileName));
    if LoadData == true && isfile(MatFilePath)
        load(MatFilePath);
    else
        [Data,H,CH] = DatRead(FilePath,'forcereading','true','datatype','ushort');
        save(MatFilePath,'Data','CH','H');
    end
    if ndims(Data) == 3
        Data = permute(Data,[4 1 2 3]);
    end
    %% Preprocess data
    SaturatedCurves = CheckSaturation(Data,coeff_DNL,OPTODE,FileName);
    if isPlot == true && LevelIsPlot >= LITE && ~isempty(SaturatedCurves)
        figure('Name',strcat(FileName,'-Saturated Curves'));
        semilogy(SaturatedCurves');
    end
    AcqTime = CH.McaTime;
    Factor = CH.McaFactor;
    NumChan = size(Data,4);
    NumRep = size(Data,1);
    ValidRepRange = round(NumRep*2/3):NumRep;
    NumValidRep = numel(ValidRepRange);
    Seq = Seq(Seq(:,1) ~= 0,:);
    FinalDummyLinesId = find(Seq(:,1) == FinalDummyLineID);
    if isempty(FinalDummyLinesId), errordlg({strcat('File: ', FileName),strcat('Add the FinalDummyLineID ',num2str(FinalDummyLineID),' to the sequence')});return; end
    Seq(FinalDummyLinesId,:) = [];
    DummyLines = find(Seq(:,1) == DummyLineID);
    if isempty(DummyLines), errordlg({strcat('File: ', FileName),strcat('Add the DummyLineID ',num2str(DummyLineID),' to the sequence')});return; end
    %     [~,~,Indexes]=unique(Seq(:,3:5),'rows','stable');
    %     DummyLines = vertcat(find(diff(Indexes) == 0));
    Seq(DummyLines,:) = [];
    Data = squeeze(Data(:,:,OPTODE,:));
    if NumRep == 1
        Data = permute(Data,[3 1 2]);
    end
    Data = permute(Data,[3 2 1]);
    Data = Data(:,:,ValidRepRange);
    DummyLines = vertcat(DummyLines,FinalDummyLinesId); %#ok<AGROW>
    Data(:,DummyLines,:) = [];
    coeff_DNL = coeff_DNL(ValidChanRange);
    NumSeqLines = size(Seq,1);
    [Delays,SortingIndexes] = unique(Seq(:,[3 4]),'rows','sorted');
    [~,Indexes]=ismember(Delays(:,1),DelayCourse(:,1));
    DelayValues = DelayCourse(Indexes,2);
    [~,Indexes]=ismember(Delays(:,2),DelayFine(:,1));
    DelayValues(:,2) = DelayFine(Indexes,2);
    DelayValues = sum(DelayValues,2);
    NumDelays = size(Delays,1);
    NumBlocks = NumSeqLines/NumDelays;
    MaskWidthCh = round(MaskWidthPs/Factor);
    
    % Riordino della sequenza
    SortingIndexes = repmat(SortingIndexes,NumBlocks,1)+...
        repelem((0:NumDelays:NumSeqLines-NumDelays)',NumDelays);
    Seq = Seq(SortingIndexes,:);
    Data = Data(:,SortingIndexes,:);
    %
    % Selezione dei delays
    if isinf(DelayRange)
        ValiDelayValues = 1:NumDelays;
        ValiDelayRange = 1:NumSeqLines;
    else
        ValiDelayValues = DelayRange;
        ValiDelayRange = repmat(DelayRange',NumBlocks,1)+repelem((0:NumDelays:NumSeqLines-NumDelays)',numel(DelayRange));
    end
    Seq = Seq(ValiDelayRange,:);
    Data = Data(:,ValiDelayRange,:);
    NumSeqLines = size(Seq,1);
    DelayValues = DelayValues(ValiDelayValues,:);
    NumDelays = size(Delays,1);
    NumBlocks = NumSeqLines/NumDelays;
    NumLasers = NumBlocks/2;
    
    isSeqT = false;
    StepBlock = 2;
    if contains(SeqFileName,'_seqt')
        NumLasers = NumBlocks;
        StepBlock = 1;
        isSeqT = true;
    end
   
    
    LaserId = unique(Seq(:,5),'stable');
    LaserId = LaserId(1:StepBlock:end);
    
    SPADonData = squeeze(Data(AreaCh,:,:));
    SPADonData = reshape(SPADonData,NumDelays,NumBlocks,NumValidRep);
    
    IntensityData = Data(IntensityCh,:,:);
    IntensityData = reshape(IntensityData,2,NumDelays,NumBlocks,NumValidRep);
    IntensityData(1,:,:,:) = IntensityData(1,:,:,:)*power(2,16);
    Intensity = squeeze(sum(IntensityData,1));
    
    TempData = squeeze((Data(TempChan,:,:)+3000)/100);
    TempData = reshape(TempData,NumDelays,NumBlocks,NumValidRep);
    
    Curves = Data(ValidChanRange,:,:);
    Curves = permute(Curves,[2 1 3]);
    Curves = reshape(Curves,NumDelays,NumBlocks,NumValidChan,NumValidRep);
    
    SPADOn = mean(SPADonData,3);
    SPADOn = SPADOn(:,1:StepBlock:NumBlocks);
    MeanTempRepData = mean(TempData,3);
    TempRepData = squeeze(mean(TempData,[1 2]));
    SignalTemp = MeanTempRepData(:,1:StepBlock:NumBlocks);
    if StepBlock == 2
        DarkTemp = MeanTempRepData(:,(1:StepBlock:NumBlocks)+1);
    else
        DarkTemp = SignalTemp;
    end
    disp(strcat(FileName,' - Average temp per laser:',32,num2str(mean(SignalTemp,1))))
    Signal = Curves(:,1:StepBlock:NumBlocks,:,:);
    if StepBlock == 2
        Dark = Curves(:,(1:StepBlock:NumBlocks)+1,:,:);
    else
        Dark = Signal;
    end
    SignalIntensity = Intensity(:,1:StepBlock:NumBlocks,:);
    if StepBlock == 2
        DarkIntensity = Intensity(:,(1:StepBlock:NumBlocks)+1,:);
    else
        DarkIntensity = SignalIntensity;
    end
    
    if CompensateCountingLoss == true
        for il = 1 : NumLasers
            [curve_corrected_sign, coeff_corr_s] = deadtime_correction_matrix (squeeze(Signal(:,il,:,:)),NumChan,ValidChanRange(1),ValidChanRange(end),AcqTime);
            [curve_corrected_dark, coeff_corr_d] = deadtime_correction_matrix (squeeze(Dark(:,il,:,:)),NumChan,ValidChanRange(1),ValidChanRange(end),AcqTime);
            Signal(:,il,:,:) = curve_corrected_sign;
            Dark(:,il,:,:) = curve_corrected_dark;
        end
    end
    
    SignalAllRep = mean(Signal,4);
    SignalAllRepCW = sum(SignalAllRep,3);
    AP_add = repmat(SignalAllRepCW,[1,1,NumValidChan])*(AP_prob/NumValidChan);
    DarkAllRep = mean(Dark,4);
    CW_all = sum(Signal,3);
    CW_mean = squeeze(mean(CW_all,4));
    CW_std = squeeze(std(CW_all,0,4));
    CW_all = squeeze(sum(Signal,3));
    
    DiffCounts = (SignalIntensity-CW_all)./SignalIntensity * 100;
    if any(DiffCounts(:)>10)
        warning(strcat('Could be saturation in file: ',FileName))
    end
    
    if isPlot == true
        if LevelIsPlot>=LITE
            figure('Name',strcat(FileName,'-Temperature vs Delays'));
            plot(1:NumDelays*NumLasers,[SignalTemp(:),DarkTemp(:)])
            legend({'Signal' 'Dark'})
            if NumRep ~= 1
                figure('Name',strcat(FileName,'-Temperature vs Repetition'));
                plot(ValidRepRange,TempRepData)
            end
        end
        
        if LevelIsPlot>=HARD
            figure('Name',strcat(FileName,'-SIGNAL. Raw curves (repetition summed)'));
            %         tiledlayout('flow')
            nsub = numSubplots(NumLasers);
            for il = 1:NumLasers
                %           nexttile;
                subplot(nsub(1),nsub(2),il);
                semilogy(ValidChanRange,squeeze(SignalAllRep(:,il,:)));
                title(strcat('LD #:',32,num2str(LaserId(il))));
                grid on;
                %             xlabel ('Channel (a.u.)');
                %             ylabel ('Counts (a.u.)');
                xlim([20 ValidChanRange(end)]);
                ylim([1 10^5])
            end
            
            figure('Name',strcat(FileName,'-DARK. Raw curves (repetition summed)'));
            %         tiledlayout('flow')
            for il = 1:NumLasers
                %             nexttile;
                subplot(nsub(1),nsub(2),il);
                semilogy(ValidChanRange,squeeze(DarkAllRep(:,il,:)));
                title(strcat('LD #:',32,num2str(LaserId(il))));
                grid on;
                %             xlabel ('Channel (a.u.)');
                %             ylabel ('Counts (a.u.)');
                xlim([20 ValidChanRange(end)]);
                ylim([1 10^5])
            end
        end
        
        if LevelIsPlot>= MEDIUM
            for il = 1:NumLasers
                figure('Name',strcat(FileName,'-Laser id: ',32,num2str(LaserId(il))));
                %             tiledlayout('flow')
                nsub = numSubplots(NumDelays);
                for id = 1:NumDelays
                    %                 nexttile;
                    subplot(nsub(1),nsub(2),id);
                    semilogy(ValidChanRange,squeeze(Signal(id,il,:,:)));
                    title(strcat('LD #:',32,num2str(LaserId(il)),'. Delay #:',32,num2str(id)));
                    grid on;
                    %                 xlabel ('Channel (a.u.)');
                    %                 ylabel ('Counts (a.u.)');
                    xlim([20 ValidChanRange(end)]);
                    ylim([1 10^5])
                end
            end
        end
    end
    
    %% Analysis
    % DNL correction
    DNL_matrix = permute(repmat(coeff_DNL,[1,NumDelays,NumLasers]),[2,3,1]);
    if DNLCorrection == true
        SignalAllRep = SignalAllRep.*DNL_matrix;
        DarkAllRep = DarkAllRep.*DNL_matrix;
    end
    
    % Remove bkg
    if RemoveBkg == true && isSeqT == false
        SignalAllRep = SignalAllRep - (DarkAllRep * BkgWeigth + AP_add) ;
    end
    
    if isPlot == true
        if LevelIsPlot >= HARD
            figure('Name',strcat(FileName,'-SIGNAL. Corrected curves (repetition summed)'));
            tiledlayout('flow')
            for il = 1:NumLasers
                ah = nexttile;
                semilogy(ValidChanRange,squeeze(SignalAllRep(:,il,:)));
                title(ah,strcat('LD:',32,num2str(LaserId(il)),'. Raw curves (repetition summed), DNL correction:',32,num2str(DNLCorrection),32,'. BKG SUB:',32,num2str(RemoveBkg)));
                grid on;
                %             xlabel ('Channel (a.u.)');
                %             ylabel ('Counts (a.u.)');
                ylim([1 10^5])
            end
            figure('Name',strcat(FileName,'-DARK. Corrected curves (repetition summed)'));
            tiledlayout('flow')
            for il = 1:NumLasers
                ah = nexttile;
                semilogy(ValidChanRange,squeeze(DarkAllRep(:,il,:)));
                title(ah,strcat('LD:',32,num2str(LaserId(il)),'. Raw curves (repetition summed), DNL correction:',32,num2str(DNLCorrection),32,'. BKG SUB:',32,num2str(RemoveBkg)));
                grid on;
                %             xlabel ('Channel (a.u.)');
                %             ylabel ('Counts (a.u.)');
                ylim([1 10^1])
            end
        end
    end
    
    subnum=numSubplots(NumLasers);
    if contains(SeqFileName,'_seqt') && isPlot && LevelIsPlot >= LITE
        figure('Name',strcat(num2str(FileName),'-Trimmered curves'));
        for il = 1:NumLasers
            ah = subplot(subnum(1),subnum(2),il);
            semilogy(ValidChanRange,squeeze(SignalAllRep(:,il,:)));
            title(ah,strcat('LD:',32,num2str(LaserId(il)),32,num2str(FileName)));
            grid on;
        end
        continue;
    end
    
    CurveOutput = zeros(NumChanToPad,NumLasers);
    FWHMCurveOutput = zeros(il,1);
    
    fh(1) = figure(20+ifn-1);
    fh(1).Name = strcat(num2str(FileName),'-Curve rescaled in time and amplitude');
    fh(1).WindowState = 'maximized';
    fh(2) = figure(50+ifn-1);
    fh(2).Name = strcat(num2str(FileName),'-Curve masked');
    fh(2).WindowState = 'maximized';
    fh(3) = figure(70+ifn-1);
    fh(3).Name = strcat(num2str(FileName),'-Curve reconstructed');
    fh(3).WindowState = 'maximized';
    
    for il = 1:NumLasers
        SignalSingleLaser = squeeze(SignalAllRep(:,il,:));
        SpadSingleLaser = SPADOn(:,il);
        SingleLaserCurve = padarray(SignalSingleLaser,[0,NumChanToPad-NumValidChan],0,'post');
        
        AreaOn = SPAD_meas(SpadSingleLaser);
        FWHM_ch = arrayfun(@(id) CalcWidth(SignalSingleLaser(id,:),0.5),1:NumDelays);
        %         FWHM_SingleLasers(il,:) = FWHM_ch;
        GatesShifts = round(DelayValues/Factor);
        SingleLaserCurveShifted = arrayfun(@(id) circshift(squeeze(SingleLaserCurve(id,:)),GatesShifts(id)),1:NumDelays,'UniformOutput',false)';
        SingleLaserCurveShifted = cell2mat(SingleLaserCurveShifted);
        
        CoeffNorm = AreaOn./AreaOn(end);
        CoeffNorm = repmat(CoeffNorm, [1,NumChanToPad]);
        SingleLaserCurveShiftedRescaled = SingleLaserCurveShifted./CoeffNorm;
        
        set(0, 'CurrentFigure', fh(1))
        subplot(subnum(1),subnum(2),il);
        semilogy(1:NumChanToPad,SingleLaserCurveShiftedRescaled);
        title(strcat('LD:',32,num2str(LaserId(il)),32,num2str(FileName),32,'.Curve rescaled in time and amplitude'));
        grid on;
        
        % masking
        ReconCurve = ones(NumDelays,NumChanToPad) * NaN;
        
        if AutoMask == true
            MaxVal=max(SingleLaserCurveShiftedRescaled(1,:));
            MaskStart_ch(il) = find(SingleLaserCurveShiftedRescaled(1,:)>MaxVal*PercMaxAutoMask,1,'last');
        end
        
        for id = 1:NumDelays
            if id == 1
                open_ch = 1;
                close_ch = (MaskStart_ch(il)+MaskWidthCh);
            else
                open_ch = MaskStart_ch(il) + round((DelayValues(id) - DelayValues(1))/Factor) + 1;
                close_ch  =  open_ch + MaskWidthCh;
                close_ch = min(close_ch,NumChanToPad);
            end
            ReconCurve(id,open_ch:close_ch)= SingleLaserCurveShiftedRescaled(id,open_ch:close_ch);
        end
        
        CurveOutput(:,il) = mean(ReconCurve,1,'omitnan');
        CurveOutput(isnan(CurveOutput(:,il)),il) = 0;
        FWHMCurveOutput(il) = CalcWidth(CurveOutput(:,il),0.5)*Factor;
        
        set(0, 'CurrentFigure', fh(2))
        subplot(subnum(1),subnum(2),il);
        semilogy(1:NumChanToPad,ReconCurve,'linewidth',2);
        grid on;
        title(strcat('LD:',32,num2str(LaserId(il)),'. Curve masked'));
        %         xlabel('Channel (a.u.)');
        %         ylabel('Counts (a.u.)');
        xlim([30 130]); ylim([1 10^8]);
        
        set(0, 'CurrentFigure', fh(3))
        subplot(subnum(1),subnum(2),il);
        semilogy(1:NumChanToPad,(squeeze(CurveOutput(:,il))./max(CurveOutput(:,il))),col(1,1),'linewidth',2);
        grid on;
        title(strcat('LD:',32,num2str(LaserId(il)),'.FWHM:',num2str(FWHMCurveOutput(il),'%.0f'),' ps'));
        %         xlabel('Channel (a.u.)');
        %         ylabel('Counts (a.u.)');
        xlim([30 130]); ylim([10^(-5) 1]);
    end
    
    %% plot counts
    if LevelIsPlot >= HARD
        figure('Name',strcat(FileName,'. DelaysVsIntensity'))
        errorbar(CW_mean,CW_std);
        xlabel('Delays');
        ylabel('Intensity')
        grid on;
    end
    
    if IsWriteCurveToFile == true
        load(strcat(dirdata,'templatesub.mat'));
        FilePathDat = strcat(dirdata,strcat(FileName,'_recon'),'.dat');
        CH.LoopNum(1:2) = [1 1];
        CH.McaChannNum = NumChanToPad;
        H = CompileHeader(CH);
        fid = fopen(FilePathDat,'wb');
        fwrite(fid,H,'uint8');
        for il = 1:NumLasers
            templatesub.Fiber = il;
            templatesub.Det = il;
            SUBH = CompileSubHeader(templatesub);
            fwrite(fid,SUBH,'uint8');
            fwrite(fid,CurveOutput(:,il),'uint32');
        end
        fclose(fid);
    end
end
warning('on','MATLAB:Axes:NegativeDataInLogAxis');