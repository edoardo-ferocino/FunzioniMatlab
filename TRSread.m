function [Sett,varargout] = TRSread(SettingsName)
%TRSread(SettingsName) creates variables in workspace from .TRS file.
%SettingsName is either a path or the filename. Default path is %\Settings
%[A,B,C,D,E] = TRSread(SettingsName)
%A: MostUsed B: Type C: Type D: Numeric E: String

if ~isempty(fileparts(SettingsName))
    [Path,SettingsName,~]=fileparts(SettingsName);
else
    Path = fullfile(cd(cd('..')),'Settings');
end
SourcePath =fullfile(Path,SettingsName);
copyfile([SourcePath '.TRS'],[SourcePath '_temp.TRS']);
DestinationPath = fullfile(cd,[SettingsName '.txt']);
movefile([SourcePath '_temp.TRS'],DestinationPath);
data=readtable(DestinationPath,'Delimiter','tab','ReadVariableNames',0);
delete(DestinationPath);
%Panel = double(data.Var1(:));
Type = double(data.Var2(:));
Label = data.Var3(:);
Values = data.Var4(:);
Numeric = zeros(length(Values),1);
String = cell(length(Values),1);

MostUsed{1,1} = 'ChannNum'; MostUsed{2,1} = 'LoopNum_1'; MostUsed{3,1} = 'LoopNum_2'; MostUsed{4,1} = 'LoopNum_3'; MostUsed{5,1} = 'LoopNum_4'; MostUsed{6,1} = 'LoopNum_5';
MostUsed{7,1} = 'SpcFactor'; MostUsed{8,1} = 'SpcRoutingBits'; MostUsed{9,1} = 'SpcType'; MostUsed{10,1} = 'SpcTimeM'; MostUsed{11,1} = 'ScFirstBin'; MostUsed{12,1} = 'ScLastBin';
MostUsed{13,1} = 'SpcRepRate'; MostUsed{14,1} = 'SpcRefolding'; MostUsed{15,1} = 'SelRepNum';
for iss=1:10
        MostUsed{16+4*(iss-1),1} = strcat('TrimStep_',num2str(iss));
        MostUsed{17+4*(iss-1),1} = strcat('TrimGoal_',num2str(iss));
        MostUsed{18+4*(iss-1),1} = strcat('TrimRegion_',num2str(iss));
        MostUsed{19+4*(iss-1),1} = strcat('TrimTime_',num2str(iss));
end
NumArgOut = nargout;
NumEntries = length(Values);
NumMostUsed = length(MostUsed);
for it=1:NumEntries
    
    if Type(it)==3
        String{it} = Values{it};
        Val = String{it};
    else
        Numeric(it)=str2double(Values{it});
        Val = Numeric(it);
    end
    
    Sett.(Label{it}) = Val;
    for im=1:NumMostUsed
        if(strcmp(MostUsed{im,1},Label{it}))
            MostUsed{im,2} = str2double(Values{it});
            if(strcmp(MostUsed{im,1},'SpcRoutingBits'))
                MostUsed{im,1} = 'NumChannels';
                MostUsed{im,2} = 2^str2double(Values{it});
            end
            if(strcmp(MostUsed{im,1},'ChannNum'))
                MostUsed{im,1} = 'NumBin';
            end
            if(strcmp(MostUsed{im,1},'ScFirstBin'))
                MostUsed{im,1} = 'FirstBin';
            end
            if(strcmp(MostUsed{im,1},'ScLastBin'))
                MostUsed{im,1} = 'LastBin';
            end
            if(strcmp(MostUsed{im,1},'SpcTimeM'))
                MostUsed{im,1} = 'AcqTime';
            end
            if(strcmp(MostUsed{im,1},'SpcFactor'))
                MostUsed{im,1} = 'Factor';
            end
            if(strcmp(MostUsed{im,1},'SpcRepRate'))
                MostUsed{im,1} = 'RepRate';
            end
            if(strcmp(MostUsed{im,1},'SpcType'))
                switch str2double(Values{it})
                    case 10
                        MostUsed{im,2} = 'TDC';
                    case 6
                        MostUsed{im,2} = 'Becker';
                end
            end
            if(strcmp(MostUsed{im,1},'SpcRefolding'))
                switch str2double(Values{it})
                    case -1
                        MostUsed{im,2} = 'RawData';
                    case 0
                        MostUsed{im,2} = 'Linearized';
                    case 1
                        MostUsed{im,2} = 'Antonio';
                    case 2
                        MostUsed{im,2} = 'Alberto';
                    otherwise
                        MostUsed{im,2} = 'HardwareRefold';
                end
                MostUsed{im,1} = 'RefoldType';
            end
            Sett.(MostUsed{im,1}) = MostUsed{im,2};
        end
    end
    if strcmp(Label{it},'RoiPage_1_1')
        Roi = zeros(96,1);
        for im=1:96
            Roi(im) = str2double(Values{it+im-1});
        end
        Roi=vec2mat(Roi',3);
        Sett.Roi = Roi;
    end
end

switch NumArgOut-1
    case 1
        varargout{1} = MostUsed;
    case 2
        varargout{1} = MostUsed;
        varargout{2} = Label;
    case 3
        varargout{1} = MostUsed;
        varargout{2} = Label;
        varargout{3} = Type;
    case 4
        varargout{1} = MostUsed;
        varargout{2} = Label;
        varargout{3} = Type;
        varargout{4} = Numeric;
    case 5
        varargout{1} = MostUsed;
        varargout{2} = Label;
        varargout{3} = Type;
        varargout{4} = Numeric;
        varargout{5} = String;
end
itr=1; jump=0;
for im=1:NumMostUsed
    if((~isempty(strfind(MostUsed{im,1},'TrimStep_')))&&MostUsed{im,2}==-1)
        itr=itr+1;
        jump=4;
    else
        if(jump==0)
            %assignin('caller',MostUsed{im,1},MostUsed{im,2});
        else 
            jump=jump-1;
        end
    end
end
TrimRoi=zeros(10-(itr-1),3);
itr=1;
for im=1:NumMostUsed
    if((~isempty(strfind(MostUsed{im,1},'TrimStep_')))&&MostUsed{im,2}~=-1)
        TrimRoi(itr,:)=Roi(MostUsed{im+2,2}+1,:);
        itr=itr+1;
    end
end

%assignin('caller','Roi',Roi);
%assignin('caller','TrimRoi',TrimRoi);
end

