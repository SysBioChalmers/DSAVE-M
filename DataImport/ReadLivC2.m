function s = ReadLivC2(dir)
% ReadLivC2
%   Reads liver cancer data from file. 
% Input:
%   dir                 The directory where the files are
%
% Output: 
%   ds                  SCDataset containing the cells from the tumors
%
% Usage: ds = ReadLivC2('../../ImportableData/LiverT2')
%
% Johan Gustafsson, 2019-05-20
%

%for testing:
%dir = '../../ImportableData/LiverT2';

s = Read10xMatrix(dir, 'GSE140228_UMI_counts_Droplet.mtx', ...
                  'GSE140228_UMI_counts_Droplet_genes.tsv', ...
                  'GSE140228_UMI_counts_Droplet_barcodes.tsv');



%read and process metadata
t = readtable(strcat(dir, '/GSE140228_UMI_counts_Droplet_cellinfo.tsv'), 'ReadVariableNames',true, 'ReadRowNames', false, 'Delimiter', '\t', 'FileType', 'text');

%check if the ids come in the same order:
%ids1 = s.cellIds;
%ids2 = table2cell(t(:,1));
%sum(strcmp(ids1,ids2.')) %this gives the same number as the number of
%cells, so they are in the same order!

%patient = table2cell(t(:,2));
%tissue = table2cell(t(:,3));
celltypeSub = table2cell(t(:,4));
celltypeSub3first = extractBetween(celltypeSub,1,3);

sample = table2cell(t(:,7)); % a combination of patient and tissue it seems

s.sampleIds = sample.';
s.extraCellInfo = celltypeSub.';

numCells = size(sample,1);

temp = zeros(1,numCells);
%Get the cell type for all cells.
for i = 1:numCells
    cts = celltypeSub{i};
    cts3 = celltypeSub3first{i};

    if strcmp(cts,'CD4-c6-FOXP3')
        temp(1,i) = Celltype.TCellReg;
    elseif strcmp(cts3,'CD4')
        temp(1,i) = Celltype.TCellCD4Pos;
    elseif strcmp(cts3,'CD8')
        temp(1,i) = Celltype.TCellCD8Pos;
    elseif strcmp(cts3,'NK-')
        temp(1,i) = Celltype.NKCell;
    elseif strcmp(cts3,'Mon')
        temp(1,i) = Celltype.Monocyte;
    elseif strcmp(cts3,'DC-')
        temp(1,i) = Celltype.Dendritic;
    elseif strcmp(cts3,'Mφ') 
        temp(1,i) = Celltype.Macrophage;
    elseif strcmp(cts3,'Mas') 
        temp(1,i) = Celltype.Mast;
    elseif strcmp(cts3,'ILC') 
        temp(1,i) = Celltype.InnateLymphoid;
    elseif strcmp(cts3,'B c') 
        temp(1,i) = Celltype.BCell;
    elseif strcmp(cts3,'pla') 
        temp(1,i) = Celltype.BCell;
    elseif strcmp(cts,'Lymphoid-B') 
        temp(1,i) = Celltype.BCell;
    elseif strcmp(cts,'Lymphoid-B-Plasma') 
        temp(1,i) = Celltype.BCell;
    else
        error(strcat('ReadLivC2: Invalid cell type:', cts));
        temp(1,i) = Celltype.Unknown;
    end
end
%now assign the cell types to the right rows, needs to be mapped by cell id!
s.cellType = temp;

end
