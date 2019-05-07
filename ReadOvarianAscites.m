%Reads ovarian ascites datasetfrom file
function ds = ReadOvarianAscites(pathMatrix, pathSamples, pathGenes, pathClassification)

ds = SCDataset;
ds.name = 'ovasc';
ds.data = dlmread(pathMatrix, ',');

ts = readtable(pathSamples, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t');
ds.sampleIds = table2cell(ts(:, 1));
ds.sampleIds = ds.sampleIds.';

tg = readtable(pathGenes, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t');
ds.genes = table2cell(tg(:, 1));

%get cell classifications
a = load(pathClassification);

cts = a.predicted_celltype(a.id_ascites);
ds.cellType = arrayfun(@ImpId2CellTypeId, cts);% a tsne plot confirms that they come in the right order

ds = ds.fillEmpties();

function ret = ImpId2CellTypeId(in)
    switch in
        case 0
            ret = Celltype.Unknown;
        case 1
            ret = Celltype.TCell;
        case 2
            ret = Celltype.TCellCD4Pos;
        case 3
            ret = Celltype.TCellCD8Pos;
        case 4
            ret = Celltype.TCellReg;
        case 5
            ret = Celltype.BCell;
        case 6
            ret = Celltype.MacrophageOrMonocyte;
        case 7
            ret = Celltype.Dendritic;
        case 8
            ret = Celltype.NKCell;
        case 9
            ret = Celltype.Endothelial;
        case 10
            ret = Celltype.Fibroblast;
        case 11
            ret = Celltype.OvarianCarcinoma;
        case 12
            ret = Celltype.Melanoma;
        otherwise
            ret = Celltype.Unknown;
    end
end

end
