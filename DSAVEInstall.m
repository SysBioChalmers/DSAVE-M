%Run DSAVEInstall.Install() to install (will set up the path in matlab)
%Run DSAVEInstall.Uninstall() to clear the path from matlab
classdef DSAVEInstall
    methods (Static)
        function Install
            mainDir = fileparts(which(mfilename));
            %Add paths of all directories
            addpath(mainDir);
            addpath(strcat(mainDir,'/DataImport'));
            addpath(strcat(mainDir,'/DSAVE'));
            addpath(strcat(mainDir,'/FigureGeneration'));
            addpath(strcat(mainDir,'/GeneralFunctions'));
            addpath(strcat(mainDir,'/ProgressBar'));
            addpath(strcat(mainDir,'/Tests'));
            savepath;
        end
        function Uninstall
            mainDir = fileparts(which(mfilename));
            %Remove paths of all directories
            rmpath(mainDir);
            rmpath(strcat(mainDir,'/DataImport'));
            rmpath(strcat(mainDir,'/DSAVE'));
            rmpath(strcat(mainDir,'/FigureGeneration'));
            rmpath(strcat(mainDir,'/GeneralFunctions'));
            rmpath(strcat(mainDir,'/ProgressBar'));
            rmpath(strcat(mainDir,'/Tests'));
            savepath;
        end
    end
end
