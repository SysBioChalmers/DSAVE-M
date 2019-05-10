%Reads PBMCB10000 data from file
%do not add a slash at the end of directoryPath
function ds = ReadPBMCB10000(directoryPath)
%function ds = ReadPBMCB10000()
%directoryPath = 'C:/Work/MatlabCode/components/SCLib/ImportableData/PBMC10000BCells/filtered_matrices_mex/hg19';

ds = Read10xMatrix(directoryPath);
ds.cellType(1,:) = Celltype.BCell;
ds.name = 'pbmc b 10000';

end
