function [FileName,PathName,FilterId]=uiputfilecustom(varargin)
AppDataTempPath = getapplicationdatadir(fullfile('\Temp'),false,true);
TempPath = fullfile(AppDataTempPath,'temppath.mat');
if exist(TempPath,'file')    % Load the Sim.mat file
    load(TempPath,'PathName');
    switch numel(varargin)
        case 0
            Filter = '*'; Title = 'Choose file'; 
        case 1
            Filter = varargin{1}; Title = 'Choose file'; 
        case 2
            Filter = varargin{1}; Title = varargin{2}; PathName = pwd;
    end
    [FileName,PathName,FilterId] = uiputfile(Filter,Title,PathName);
    
else
    [FileName,PathName,FilterId] = uiputfile(varargin{:});
    
end
if PathName == 0, return, end
VarArgin = varargin;
if isempty(VarArgin)
    save(TempPath,'PathName');
else
    save(TempPath,'PathName','VarArgin');
end
end