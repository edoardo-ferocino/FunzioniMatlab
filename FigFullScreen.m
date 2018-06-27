function [ fig_handle ] = FigFullScreen( varargin )
Prop=java.awt.Toolkit.getDefaultToolkit().getScreenSize();
scrsz = get(0,'ScreenSize');
pos = 1000*[0.0010    0.0410    1.5360    0.7488];
pos = [1 1 1536 788.800000000000];
% if nargin==1
%     fig_handle = figure('NumberTitle','off','Name',varargin{1},'Position',[0 0 scrsz(3) scrsz(4)]);
% else 
%     fig_handle = figure('Position',[0 0 scrsz(3) scrsz(4)]);
% end
if nargin==1
    fig_handle = figure('NumberTitle','off','Name',varargin{1},'Position',pos);
else 
    fig_handle = figure('Position',pos);
end
end

