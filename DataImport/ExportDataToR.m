%Exports the datasets used to build the standard template to R

ovm = SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte);
bc2t = SCDep.scd_bc2.cellSubset(SCDep.scd_bc2.cellType == Celltype.TCellCD4Pos | SCDep.scd_bc2.cellType == Celltype.TCellCD8Pos | SCDep.scd_bc2.cellType == Celltype.TCellReg);
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
b10000 = SCDep.scd_pbmcb10000;
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

% datasets = {ovm,bc2t_bc4tumor, b10000, scd_GSE112845_cd8};

mainDir = fileparts(which(mfilename));
cd(mainDir);
if(~exist('../../ExportedData','dir'))
    mkdir('../../ExportedData');
end


ovm.saveDataTable('../../ExportedData/ovm.txt');

numIds = 1:size(bc2t_bc4tumor.cellIds,2);
bc2t_bc4tumor.cellIds = arrayfun(@(x) ['c' num2str(x)],numIds,'UniformOutput',false);
bc2t_bc4tumor.saveDataTable('../../ExportedData/bc2t_bc4tumor.txt');

[~, ggg] = unique(b10000.genes);
b10000 = b10000.geneSubset(b10000.genes(ggg,:));%just throw away gene duplicates;
numIds = 1:size(b10000.cellIds,2);
b10000.cellIds = arrayfun(@(x) ['c' num2str(x)],numIds,'UniformOutput',false);
b10000.saveDataTable('../../ExportedData/b10k.txt');

[~, ggg] = unique(scd_GSE112845_cd8.genes);
scd_GSE112845_cd8 = scd_GSE112845_cd8.geneSubset(scd_GSE112845_cd8.genes(ggg,:));%just throw away gene duplicates;
numIds = 1:size(scd_GSE112845_cd8.cellIds,2);
scd_GSE112845_cd8.cellIds = arrayfun(@(x) ['c' num2str(x)],numIds,'UniformOutput',false);
scd_GSE112845_cd8.saveDataTable('../../ExportedData/GSE112845_cd8.txt');
