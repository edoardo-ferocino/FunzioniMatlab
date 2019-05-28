function [Width, varargout]=CalcWidth(Data,Level)
%[Width,Baricentre,MaxPos,MaxVal]=CalcWidth(Data,Level);
%Level is percentage from base of the pulse

NumChan=length(Data);
[Peak,IndexPeak]=max(Data);
DataNorm=Data/Peak;
right=find(DataNorm(IndexPeak:end)<Level,1,'first')+(IndexPeak-1);
if isempty(right), right=NumChan; end
left=find(DataNorm(1:IndexPeak)<Level,1,'last');
if isempty(left), left=1; end
extra=(Level-DataNorm(left))/(DataNorm(left+1)-DataNorm(left))+(Level-DataNorm(right))/(DataNorm(right-1)-DataNorm(right));
Width=(right-left)-extra;
NumArgOut = nargout-1;
if NumArgOut>=1
        right=find(DataNorm(IndexPeak:end)<0.05,1,'first')+(IndexPeak-1);
        if isempty(right), right=NumChan; end
        left=find(DataNorm(1:IndexPeak)<0.05,1,'last');
        if isempty(left), left=1; end
        somma =0;
        for ib=left:right
            somma=somma +Data(ib)*ib;
        end
        varargout{1}= somma/sum(Data(left:right));
end
if nargout >=2
        [A,varargout{2}] = max(Data); 
end
if nargout >=3
    varargout{3} = A;
end
end
