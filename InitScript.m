TempPath = getapplicationdatadir(fullfile('\Temp'),false,true);
TempFile = fullfile(TempPath,'tempmatlab.mat');
if exist(TempFile,'file')
    delete(TempFile)
end
clearvars
close all
clc