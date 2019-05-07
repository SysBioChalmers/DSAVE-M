%Reads PBMCTCDMem10000 data from file
%This dataset comes from the publication: Zheng et al, “Massively parallel digital transcriptional profiling of single cells”
%do not add a slash at the end of directoryPath
function ds = ReadPBMCTCD4Mem10000(directoryPath)

ds = Read10xMatrix(directoryPath);
ds.cellType(1,:) = Celltype.TCellCD4Memory;
ds.name = 'pbmc T CD4Mem 10000';

end
