function MinimizeFFS(FH)
W='MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame';
warning('off',W);
jFrame = get(handle(FH),'JavaFrame');
if jFrame.isMaximized
    jFrame.setMaximized(false)
end
warning('on',W);
end