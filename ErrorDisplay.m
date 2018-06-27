function ErrorDisplay( ME )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Ns = numel(ME.stack);
if Ns>=2
    errordlg({ME.message,'in function: ',ME.stack(1).name,[' at line: ' num2str(ME.stack(1).line)],ME.stack(2).name,[' at line: ' num2str(ME.stack(2).line)]},'Error','nonmodal');
else
    errordlg({ME.message,'in function: ',ME.stack(1).name,[' at line: ' num2str(ME.stack(1).line)]},'Error','nonmodal');
end
end

