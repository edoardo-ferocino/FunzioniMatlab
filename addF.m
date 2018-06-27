function F = addF(iF,Label,Unit,minV,maxV,Levels,Treshold,F)
% addF(Label,Unit,F): Add one entry in the Figure of Merit (FOM) space
%
%   iP = index of the F entry  
%   Label = label used for rapresentation
%   Unit = measurement unit of the value
%   F = F structure

F(iF).Unit=Unit;
if(isempty(F(iF).Label)==0), 
    disp('Error in addF: Duplicated iF'); 
end
F(iF).Label=Label;
F(iF).Treshold=Treshold;
F(iF).minV=minV;
F(iF).maxV=maxV;
F(iF).Levels=Levels;
nF=numel(F);
if(nF~=iF), 
    disp('Error in addF: Mismatch beetween nF and iF'); 
end
for i=1:nF
    if(isempty(F(i).Label)==1), 
        disp('Error in addF: Empty Elements'); 
    end
end
end

