function [lls, genes, geneLls] = DSAVEGetSingleCellDivergence(ds, minUMIsPerCell, progrBarCtxt, TPMLowerBound, iterations)
% DSAVEGetSingleCellDivergence
%   Calculates the DSAVE cell-wise divergence metric.
% Input:
%   ds              The dataset (cell population) to be investigated
%   minUMIsPerCell  (optional) All cells are downsampled to this number.
%                   Cells that have less UMIs will have a NaN divergence.
%                   Defaults to 200 TPM/CPM.
%   progrBarCtxt    (optional) Progress bar context.
%   TPMLowerBound   (optional) Genes below this threshold are not 
%                   investigated. Note that use of this parameter reduces
%                   the number of UMIs per cell. Defaults to 0.
%   iterations      (optional) Number of times the process will be repeated.
%                   The result is the median of several runs; there is
%                   stochasticity here since cells are randomly
%                   down-sampled before the calculation is made.
%
% Output:
%   lls             Log-likelihood for getting the observed counts'
%                   distribution for each cell when sampling counts from 
%                   the mean dataset gene expression. The values are 
%                   negative, and the lower (i.e. more negative) the value, 
%                   the more divergent the cell is.
%
% Usage: lls = DSAVEGetSingleCellDivergence(ds);
%
% Johan Gustafsson, 2019-05-21
%
if nargin < 2
    minUMIsPerCell = 200;
end
if nargin < 3
    progrBarCtxt = [];
end
if nargin < 4
    TPMLowerBound = 0;
end
if nargin < 5
    iterations = 15;
end

%remove all lowly expressed genes
me = TPM(mean(ds.data,2));
ds = ds.geneSubset(me >= TPMLowerBound);

genes = ds.genes;

numCells = size(ds.data,2);
numGenes = size(ds.data,1);

%Figure out what to downsample to. If we have cells with fewer UMIs than
%minUMIsPerCell, those can not be evaluated
origUMIs = sum(ds.data,1);
UMIsPerCell = origUMIs;
targetUMIsPerCell = max(min(origUMIs,[],2),minUMIsPerCell);

%Get mean expression and create probabilities
expr = TPM(ds.data);
meanRefExpr = mean(expr,2);
%meanRefExpr(isnan(meanRefExpr)) = 0;
prob = meanRefExpr./sum(meanRefExpr,1);
%unsparse...
prob = full(prob);

allLls = zeros(iterations,numCells);

allGeneLls = zeros(iterations,numGenes,numCells);

progbar = ProgrBar(['GetSingleCellDivergence ''' ds.name ''''], progrBarCtxt);

logFacs = GetLogFactorials(targetUMIsPerCell);

%run several times since there is randomness in downsampling
for it = 1:iterations
    %downsample to the right number of UMIs per cell
    dsd = ds;
    toRemUMIs = UMIsPerCell - targetUMIsPerCell;
    toRemUMIs(toRemUMIs < 0) = 0;
    for i = 1:numCells
        %first, randomly select a number of indexes from the total UMIs to
        %remove
        indexesToRem = randperm(origUMIs(1,i), toRemUMIs(1,i));
        %Then create index edges for each gene, so if a gene has 5 UMIs the
        %edges to the left and right will differ 5.
        cellData = dsd.data(:,i);
        edges = [0;cumsum(cellData)+0.1];%we need to add 0.1 to the edges since histcounts checks if edge(k) <= index < edge(k+1). Otherwise, index 1 would not end up between 0 and 1
        %Now get the number of index hits you get within the edge range for
        %each gene
        [subtr,~] = histcounts(indexesToRem,edges);
        dsd.data(:,i) = dsd.data(:,i) - subtr.';
    end

    %loop through all cells and calculate log likelihood
    for i = 1:numCells
        if UMIsPerCell(1,i) < targetUMIsPerCell
            allLls(it,i) = NaN;
            allGeneLls(it,:,i) = NaN;
        else
            %Use multinomial distribution pdf
            allLls(it,i) = LogMultinomialPDF(dsd.data(:,i), prob);
            %Also fill in gene binomial pdf
            allGeneLls(it,:,i) = LogBinomialPDF(dsd.data(:,i), prob, logFacs);
        end
    end
    progbar.Progress(it/iterations);
end
%take the median of all runs
lls = median(allLls,1);

%replace all nan in the gene-wise with 0. So, nan means that the gene is not
%expressed, and then we are certain of the outcome -> PDF = 1
%->log(PDF) = 0
allGeneLls(:,isnan(prob),:) = 0;

geneLls = reshape(median(allGeneLls,1), numGenes, numCells);
progbar.Done();

end
