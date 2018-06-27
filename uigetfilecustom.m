function [FileName,PathName,FilterIndex]=uigetfilecustom(varargin)
AppDataTempPath = getapplicationdatadir(fullfile('\Temp'),false,true);
TempPath = fullfile(AppDataTempPath,'temppath.mat');
if exist(TempPath,'file')    % Load the Sim.mat file
    load(TempPath,'PathName');
    switch numel(varargin)
        case 0
            Type = {'*.m;*.fig;*.mat;*.mdl', 'All MATLAB Files'}; Title = 'Select File to Open';
        case 1
            Type = varargin{1}; Title = 'Select File to Open';
        case 2
            Type = varargin{1}; Title = varargin{2};
    end
    [FileName,PathName,FilterIndex] = uigetfile(Type,Title,PathName);
    
else
    [FileName,PathName,FilterIndex] = uigetfile(varargin{:});
    
end
if FilterIndex == 0, return, end
VarArgin = varargin;
if isempty(VarArgin)
    save(TempPath,'PathName');
else
    save(TempPath,'PathName','VarArgin');
end
OutputPath= [PathName,FileName];
end