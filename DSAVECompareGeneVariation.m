%This function is more experimental and has not been fully evaluated. The
%idea is that you can align two cell populations and compare the variation
%for each gene between them, such that to check if the variation of some
%genes have changed between conditions.
function [logCVDifference1, logCVDifference2, genes] = DSAVECompareGeneVariation(ds1, ds2, lb, maxNumCells)
if nargin < 4
    maxNumCells = 10000;
end
if nargin < 3
    lb = 10;
end
%First get rid of all genes that are not somewhat highly expressed in both
%samples. No point in looking at the other genes, it will almost always just be ones and
%zeros anyway
[ds1,ds2] = SynchronizeGenes(ds1, ds2, true);

ds1tpm = TPM(mean(ds1.data, 2));
ds2tpm = TPM(mean(ds2.data, 2));

sel = ds1tpm >= lb & ds2tpm >= lb;

ds1 = ds1.geneSubset(sel);
ds2 = ds2.geneSubset(sel);

numCells = min(min(size(ds1.data,2), size(ds2.data,2)), maxNumCells);
d1 = ds1.randSample(numCells);
d2 = ds2.randSample(numCells);


%Then downsample the cells so they have the same distribution of counts
%over the genes, also including the same number of cells
d1 = DownSampleCells(d1,d2);
d2 = DownSampleCells(d2,d1);
 

%then downsample the genes to the same expression level in both datasets
d1 = DownSampleGenes(d1,d2);
d2 = DownSampleGenes(d2,d1);

s = GenSampDs(d1);


logCV1 = GetPoolData(d1);
logCV2 = GetPoolData(d2);
logCVs = GetPoolData(s);

logCVDifference1 = logCV1 - logCVs;
logCVDifference2 = logCV2 - logCVs;
genes = ds1.genes;

function logCV = GetPoolData(ds1)

    avgRefExpr = mean(ds1.data,2);
    numRepetitions = 500;
    numGenes = size(ds1.data,1);
    numCells = size(ds1.data,2);
    poolData = zeros(numGenes, numRepetitions);
    for repetition = 1:numRepetitions
        as = datasample(ds1.data, numCells, 2);%bootstrap
        poolData(:,repetition) = TPM(sum(as,2));%take the mean over the cells -> a vertical vector with the mean value for each gene
    end
    %calc variances
    variances = var(poolData, 0, 2);
    sd = sqrt(variances);
    cv_ = sd ./ (avgRefExpr + 0.05);%Coefficient of Variation = std/mean. Adding 0.05, a neglectably small number, to handle too lowly expressed genes
    logCV = log2(cv_ + 1);%the + 1 says that no variance -> 0 value
end


function ds = DownSampleCells(d1,d2)
    %create empty dataset
    ds = d1;
    numGenes = size(ds.genes,1);
    numCells = size(ds.data,2);
    ds.data = zeros(numGenes,numCells);%sparsify in the end

    origUMIs = full(sum(d1.data,1));%sum of UMIs for each cell
    matchUMIsOrig = full(sum(d2.data,1));%sum of UMIs for each cell
    matchUMIs = matchUMIsOrig;%same size, don't randomize anything

    %sort both UMI vectors and try to downsample to match the match set
    %as good as possible
    [~,iOrig] = sort(origUMIs);
    [~,iMatchUMIs] = sort(matchUMIs);
    idealUMIs = zeros(1,numCells);
    idealUMIs(iOrig) = matchUMIs(iMatchUMIs);
    %so, only downsample the cells that are higher, leave
    %the rest. Those will be fixed when d1 and d2 switch places in the
    %subsequent call to this function
    newUMIs = min(origUMIs, idealUMIs);
    
    toRemUMIs = origUMIs - newUMIs;
    
    %Loop through all cells
    for i = 1:numCells
        %first, randomly select a number of indexes from the total UMIs to
        %remove
        indexesToRem = randperm(origUMIs(1,i), toRemUMIs(1,i));
        %Then create index edges for each gene, so if a gene has 5 UMIs the
        %edges to the left and right will differ 5.
        cellData = d1.data(:,i);
        edges = [0;cumsum(cellData)+0.1];%we need to add 0.1 to the edges since histcounts checks if edge(k) <= index < edge(k+1). Otherwise, index 1 would not end up between 0 and 1
        %Now get the number of index hits you get within the edge range for
        %each gene
        [subtr,~] = histcounts(indexesToRem,edges);
        %size(inDs.data(:,i))
        %size(subtr)
        ds.data(:,i) = d1.data(:,i) - subtr.';
    end

    %now sparsify, will be slower to work with a sparse matrix in the loop
    %ds.data = sparse(ds.data);
end
 
%will make sure each gene is expressed equally
function ds = DownSampleGenes(d1,d2)
    %create empty dataset
    ds = d1;
    numGenes = size(ds.genes,1);
    numCells = size(ds.data,2);
    ds.data = zeros(numGenes,numCells);%sparsify in the end

    origUMIs = full(sum(d1.data,2));%sum of UMIs for each gene
    matchUMIsOrig = full(sum(d2.data,2));%sum of UMIs for each gene
    matchUMIs = matchUMIsOrig;%same size, don't randomize anything

    idealUMIs = matchUMIs;
    %so, only downsample the genes that are higher, leave
    %the rest. Those will be fixed when d1 and d2 switch places in the
    %subsequent call to this function
    newUMIs = min(origUMIs, idealUMIs);
    
    toRemUMIs = origUMIs - newUMIs;
    
    %Loop through all genes
    for i = 1:numGenes
        %first, randomly select a number of indexes from the total UMIs to
        %remove
        indexesToRem = randperm(origUMIs(i,1).', toRemUMIs(i,1).');
        %Then create index edges for each gene, so if a gene has 5 UMIs the
        %edges to the left and right will differ 5.
        geneData = d1.data(i,:);
        edges = [0 (cumsum(geneData)+0.1)].';%we need to add 0.1 to the edges since histcounts checks if edge(k) <= index < edge(k+1). Otherwise, index 1 would not end up between 0 and 1
        %Now get the number of index hits you get within the edge range for
        %each gene
        [subtr,~] = histcounts(indexesToRem,edges);
        %size(inDs.data(:,i))
        %size(subtr)
        ds.data(i,:) = d1.data(i,:) - subtr;
    end
end

function ds = GenSampDs(d1)
    %create empty dataset
    ds = d1;
    numCells = size(d1.data,2);
    numGenes = size(ds.genes,1);
    ds.data = zeros(numGenes,numCells);%sparsify in the end


    %generate probabilities
    tmpTempl = TPM(d1);
    meanRefExpr = mean(tmpTempl.data,2);
    %meanRefExpr(isnan(meanRefExpr)) = 0;
    prob = meanRefExpr./sum(meanRefExpr,1);
    %unsparse...
    prob = full(prob);
    %generate a vector with the sum of probabilities up to this gene (including this gene)
    probSum = prob;%just create a vector of the right size
    probSum(1,1) = 0;%not really needed, but here for clarity
    for ii = 2:numGenes
        probSum(ii,1) = probSum(ii-1,1) + prob(ii-1,1);
    end

    edges = [probSum;max(probSum(end,1),1)];%add right edge; make sure it is not smaller than the previous due to roundoff problems. 

    %create vector of number of UMIs
    UMIs = sum(d1.data,1);

    %generate cells
    for i = 1:numCells
        if UMIs(1,i) > 0 %handle the fact that cells sometimes can contain 0 UMIs, which doesn't work with the code below. The data is automatically 0 in that case.
            %ProgressBar(floor(100*i/numCells));
            %generate n random numbers between 0-1 and sort them
            r = rand([UMIs(1,i) 1]);
            ds.data(:,i) = histcounts(r,edges);
        end
    end

    ds.data = sparse(ds.data);

end




end

