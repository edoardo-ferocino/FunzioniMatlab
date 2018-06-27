function P = addP(Order,iP,Default,Range,Label,Unit,Title,Help,P)

% addP(Order,iP,Default,Range,Label,Unit,Title,P): Add one entry in the parameter space
%
%   Order = Order in the Output (0=take default, 1=x-axis, 2=yaxis, 3=rows, 4=columns
%   iP = index of the P entry  
%   Default = default value
%   Range = range of values
%   Label = label used for rapresentation
%   Unit = measurement unit of the value
%   Dim = number of elements = 1 if Order = 0, else num elements in Range
%   P = P structure

P(iP).Order=Order;
if(isempty(P(iP).Default)==0), disp('Error in addP: Duplicated iP'); end
P(iP).Default=Default;
P(iP).Range=Range;
P(iP).Label=Label;
P(iP).Unit=Unit;
P(iP).Title=Title;
P(iP).Help=Help;
nP=numel(P);
if(nP~=iP), disp('Error in addP: Mismatch beetween nP and iP'); end
for i=1:nP
    if(isempty(P(i).Default)==1), disp('Error in addP: Empty Elements'); end
end
if Order==0, P(iP).Dim=1; else P(iP).Dim=numel(Range); end
end

