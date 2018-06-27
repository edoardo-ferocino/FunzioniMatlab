function OH = getParentFigure(H,First,Last)
% if the object is a figure or figure descendent, return the
% figure. Otherwise return [].
CF = findobj('Name',H.InterfaceName);
Childrens = CF.Children;
iO = 1;
for iH = First:Last
    for iCH = 1:numel(Childrens)
        if strcmp(H.O(iH).Tag,Childrens(iCH).Tag)
            OH(iO) = findobj('Tag',H.O(iH).Tag);
            iO = iO +1;
            break;
        end
    end
%       HL(iH) = get(findobj('tag',H.O(iH).Tag,' H.LastInstace]))
                    %set(findobj('tag',H.O(iH).Tag),'Parent',findobj('tag',H.RootInterfaceName));
end

