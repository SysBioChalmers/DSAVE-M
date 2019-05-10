%Reads 10x files into a SCDataset, i.e. the files  '/matrix.mtx',
%'/genes.tsv' and '/barcodes.tsv'
function ds = Read10xMatrix(directoryPath)
%directoryPath = 'C:/Work/MatlabCode/components/SCLib/ImportableData/PBMC10000BCells/filtered_matrices_mex/hg19';
matrixPath = strcat(directoryPath,'/matrix.mtx');
genesPath = strcat(directoryPath,'/genes.tsv');
barcodesPath = strcat(directoryPath,'/barcodes.tsv');

ds = SCDataset;

%couldn't get readmtx to work, using fscanf in combination with spconvert
%instead
%first count the number of comment rows
fileID = fopen(matrixPath);
numCommentRows = 0;
while 1
    tline = fgetl(fileID);
    if (tline(1) == '%')
        numCommentRows = numCommentRows + 1;
    else
        break;
    end
end
fclose(fileID);

%now start reading the file again
fileID = fopen(matrixPath);
%get rid of the comment rows
for i = 1:numCommentRows
    fgetl(fileID);
end

%get rid of size
fgetl(fileID);
%now read
formatSpec = '%d\t%d\t%d';
sizeA = [3 Inf];
A = fscanf(fileID,formatSpec,sizeA);
A = A.';
%max(A,[],1)
ds.data = spconvert(A);

fclose(fileID);

f = readtable(genesPath, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t', 'FileType', 'text');
ds.genes = table2cell(f(:, 2));

%fill out the last few genes with zeros, they are not part of the matrix
[matrows,cols] = size(ds.data);
generows = size(ds.genes,1);
if (matrows < generows)
    ds.data = [ds.data;zeros(generows-matrows,cols)];
end

f = readtable(barcodesPath, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t', 'FileType', 'text');
ds.cellIds = table2cell(f(:, 1)).';

ds = ds.fillEmpties();

%ds = TPM(ds);

end
