function [varargout]=save_figure(FName,varargin)
NumArgIn = nargin-1;
NumFig = 0;
ifor=1;
for in=1:NumArgIn
    if ischar(varargin{in})
        format{ifor} = varargin{in};
        ifor=ifor+1;
    end
    if ishandle(varargin{in})
        fig_handle = varargin{in};
        NumFig = length(varargin{in});
    end
    
end
if NumFig == 0
    fig_handle = gcf;
    NumFig =1;
end

for ih = 1:NumFig
set(fig_handle(ih), 'Color', 'w');
set(fig_handle(ih), 'PaperPosition', [-0.5 -0.25 6 5.5]); %Position the plot further to the left and down. Extend the plot to fill entire paper.
set(fig_handle(ih), 'PaperSize', [5 5]); %Keep the same paper size
hold off;
if iscell(FName)
    Nome = FName{ih};
else, Nome = FName;
end
 saveas(fig_handle(ih),Nome,'fig')
 %export_fig(FName{ih},'-painters','-eps','-pdf','-jpg');
 switch ifor-1
     case 0
         export_fig(Nome,'-painters','-jpg',fig_handle(ih));
     case 1
         export_fig(Nome,'-painters',format{1},fig_handle(ih));
     case 2
         export_fig(Nome,'-painters',format{1},format{2},fig_handle(ih));
     case 3
         export_fig(Nome,'-painters',format{1},format{2},format{3},fig_handle(ih));
     otherwise
         
 end
end
if nargout
varargout{1} = fig_handle;
end
end
