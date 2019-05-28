classdef DSAVEInstall
% DSAVEInstall
%   Support for installing and uninstalling DSAVE
%   Run DSAVEInstall.install() to install (will set up the path in MATLAB)
%   Run DSAVEInstall.uninstall() to clear the path from MATLAB
%
% Johan Gustafsson, 2019-05-20
%
    methods (Static)
        function install
            mainDir = fileparts(which(mfilename));
            %Add paths of all directories
            addpath(mainDir);
            addpath(strcat(mainDir,'/DataImport'));
            addpath(strcat(mainDir,'/DataExport'));
            addpath(strcat(mainDir,'/DSAVE'));
            addpath(strcat(mainDir,'/FigureGeneration'));
            addpath(strcat(mainDir,'/GeneralFunctions'));
            addpath(strcat(mainDir,'/ProgressBar'));
            addpath(strcat(mainDir,'/Tests'));
            savepath;
        end
        function uninstall
            mainDir = fileparts(which(mfilename));
            %Remove paths of all directories
            rmpath(mainDir);
            rmpath(strcat(mainDir,'/DataImport'));
            rmpath(strcat(mainDir,'/DataExport'));
            rmpath(strcat(mainDir,'/DSAVE'));
            rmpath(strcat(mainDir,'/FigureGeneration'));
            rmpath(strcat(mainDir,'/GeneralFunctions'));
            rmpath(strcat(mainDir,'/ProgressBar'));
            rmpath(strcat(mainDir,'/Tests'));
            savepath;
        end
    end
end
