function [varargout]=save_figure(FName,varargin)
NumArgIn = nargin-1;
NumFig = 0;
CustomVarargin = varargin;
for in=1:NumArgIn
    if ishandle(varargin{in})
        CustomVarargin(in) = [];
        fig_handle = varargin{in};
        NumFig = numel(fig_handle);
    end
end
if NumFig == 0
    fig_handle = gcf;
    NumFig =1;
end

for ih = 1:NumFig
    set(fig_handle(ih), 'Color', 'w');
    graphics = findobj(fig_handle(ih),'type','uipanel');
    if ~isempty(graphics)
        grcols = zeros(numel(graphics),3);
        for ig = 1:numel(graphics)
           grcols(ig,:) = graphics(ig).BackgroundColor;
           graphics(ig).BackgroundColor = 'w'; 
        end
    end
    set(fig_handle(ih), 'Color', 'w');
%     set(fig_handle(ih), 'PaperPosition', [-0.5 -0.25 6 5.5]); %Position the plot further to the left and down. Extend the plot to fill entire paper.
    set(fig_handle(ih), 'PaperSize', [5 5]); %Keep the same paper size
    set(fig_handle(ih),'PaperPositionMode','auto');
    hold off;
    if iscell(FName)
        Nome = FName{ih};
    else, Nome = FName;
    end
    savefig(fig_handle(ih),Nome)
    %export_fig(FName{ih},'-painters','-eps','-pdf','-jpg');
    if any(contains(CustomVarargin,'-svg'))
        axh=findobj(fig_handle(ih),'type','axes');
        if numel(axh)==1
            axh.LooseInset = [0 0 0.01 0];
        end
        saveas(fig_handle(ih),Nome,'svg');
        CustomVarargin(contains(CustomVarargin,'-svg')) = [];
    else
        export_fig(Nome,'-painters','-png',CustomVarargin{:},fig_handle(ih));
    end
    set(fig_handle(ih),'Color','default');
    if ~isempty(graphics)
        for ig = 1:numel(graphics)
           graphics(ig).BackgroundColor = grcols(ig,:);
        end
    end
end
if nargout
    varargout{1} = fig_handle;
end
end
