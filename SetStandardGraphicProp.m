function SetStandardGraphicProp(fh)
ah=findobj(fh,'Type','axes');
for ia=1:numel(ah)
set(ah(ia),'LineWidth',2,'FontSize',26);
Lh=findobj(ah(ia),'type','line','-or','type','constantline','-or','type','errorbar');
for il=1:numel(Lh)
    Lh(il).LineWidth = 3;
    if isprop(Lh(il),'MarkerSize')
        Lh(il).MarkerSize = 12;
    end
end
end
end