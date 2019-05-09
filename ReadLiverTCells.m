%Reads liver cancer T-cells from file into an SCDataset class
%
function  ds = ReadLiverTCells(pathData, pathCellFilter)

ds = SCDataset;
ds.name = 'LCTCells';
L = importdata(pathData,'\t');
ds.data = L.data;
[m,n] = size(ds.data);
ds.cellIds = L.textdata(1, 3:end);
ds.genes = L.textdata(2:end, 2);
%filter on cell filter, i.e. remove the cells that did not pass quality
%control in the paper:
f = readtable(pathCellFilter, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t');
okCells = table2cell(f);
[ds.cellIds, ia, ib] = intersect(ds.cellIds, okCells);
ds.data = ds.data(:, ia);

%extract patient id, format is either 'xxxx-xx-yyyy' or 'xxxx-yyyy', and we
%only want the y:s
%a bit annoying, we have to extract two tokens and only keep the last
temp = regexp(ds.cellIds, '\w+-(\w+-)?(\w+)', 'tokens');
ds.sampleIds = cellfun(@(c) c{1}{2},temp, 'UniformOutput', false);
ds.cellType = repmat(Celltype.TCell,size(ds.genes,1),size(ds.cellIds,2));

ds = ds.fillEmpties();