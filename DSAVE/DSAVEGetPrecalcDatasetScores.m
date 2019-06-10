function ret = DSAVEGetPrecalcDatasetScores()
% DSAVEGetPrecalcDatasetScores
%   Gets the DSAVE BTM variation score for 9 cell populations from 
%   different datasets. The values have been precalculated using the  
%   standard template.
%
% Usage: values = DSAVEGetPrecalcDatasetScores();
%
% Johan Gustafsson, 2019-06-10
%
    ret.names = {'BC', 'OC', 'LC', 'LIVC', 'PBMC68k', 'B10k', 'CD4TMEM', 'HCA CB', 'CD8T'};
    ret.vals = [0.1378 0.0389 0.0283 0.1102 0.0352 0.0194 0.0139 0.0163 0.0321];
end
