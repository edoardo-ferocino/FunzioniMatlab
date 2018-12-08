function SelPath=uigetdircustom(varargin)
AppDataTempPath = getapplicationdatadir(fullfile('\Temp'),false,true);
TempPath = fullfile(AppDataTempPath,'temppath.mat');
if exist(TempPath,'file')    % Load the Sim.mat file
    load(TempPath,'SelPath');
    switch numel(varargin)
        case 0
            Title = 'Select destination folder';
        case 1
            Title = varargin{1};
    end
    SelPath = uigetdir(SelPath,Title);
    
else
    SelPath = uigetdir(varargin{:});
    
end
if SelPath == 0, return, end
VarArgin = varargin;
if isempty(VarArgin)
    save(TempPath,'SelPath');
else
    save(TempPath,'SelPath','VarArgin');
end
end