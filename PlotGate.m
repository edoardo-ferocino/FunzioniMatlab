close all
FH=FFS; FH.NumberTitle = 'off'; FH.Name = 'prova';
SubH=subplot1(3,3,'Min',[0.05 0.05],'Gap',[0.05 0.05]);
for inr = 1:nW
    subplot1(inr);
    SetPushButtons(Waves,SubH,inr);
    cmin = min(Waves(inr).GatedMatrix(:));
    cmax = max(Waves(inr).GatedMatrix(:));
    imagesc(Waves(inr).GatedMatrix(:,:,1),[cmin cmax]);
    colormap pink, shading interp, axis image;
    title(num2str(Waves(inr).Wavelenghts(inr)))
    colorbar('westoutside')
end
delete(SubH(inr+1:end))
function SetPushButtons(Waves,SubH,inr)
UpControl(inr) = uicontrol; DownControl(inr) = uicontrol;
UpControl(inr).String = 'Up';DownControl(inr).String = 'Down';
UpControl(inr).Units = 'normalized';DownControl(inr).Units = 'normalized';
UpControl(inr).Position(2)=SubH(inr).OuterPosition(2)+SubH(inr).OuterPosition(4)/2;
DownControl(inr).Position(2)=SubH(inr).OuterPosition(2)+SubH(inr).OuterPosition(4)/2;
UpControl(inr).Position(1)=SubH(inr).OuterPosition(1);
DownControl(inr).Position(1)=SubH(inr).OuterPosition(1);
UpControl(inr).Units = 'pixels'; UpControl(inr).Position(2) = UpControl(inr).Position(2)+20;
DownControl(inr).Units = 'pixels';
UpControl(inr).Position(1) = UpControl(inr).Position(1)+10;
DownControl(inr).Position(1) = DownControl(inr).Position(1)+10;
UpControl(inr).Tag = ['U' num2str(inr)];DownControl(inr).Tag = ['D' num2str(inr)];
UpControl(inr).UserData = 1;DownControl(inr).UserData = 1;
UpControl(inr).Callback = {@ChangeCallback,SubH,Waves};
DownControl(inr).Callback = {@ChangeCallback,SubH,Waves};
end
function ChangeCallback(src,event,SubH,Waves)
axnum = str2double(src.Tag(2));
if strcmpi(src.String,'up')
    if src.UserData>=Waves(axnum).NumGates
        newVal = src.UserData; src.UserData = newVal;
    else
        newVal = src.UserData+1; src.UserData = newVal;
    end
    Complementary = findobj('Tag',['D' num2str(axnum)]);
    Complementary.UserData = newVal;
else
    if src.UserData<=1
        newVal = src.UserData; src.UserData = newVal;
    else
        newVal = src.UserData-1; src.UserData = newVal;
    end
    Complementary = findobj('Tag',['U' num2str(axnum)]);
    Complementary.UserData = newVal;
end
cmin = min(Waves(axnum).GatedMatrix(:));
cmax = max(Waves(axnum).GatedMatrix(:));
imagesc(SubH(axnum),Waves(axnum).GatedMatrix(:,:,newVal));
colormap pink, shading interp, axis image;
cb = colorbar(SubH(axnum),'westoutside'); cb.Limits = [cmin cmax];
title(SubH(axnum),[num2str(Waves(axnum).Wavelenghts(axnum)) '.gate ' num2str(newVal) 'of' num2str(Waves(axnum).NumGates)])
end