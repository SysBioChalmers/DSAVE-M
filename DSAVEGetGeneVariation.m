function [logCVDifference, genes, SNOVariances, pVals] = DSAVEGetGeneVariation(ds, lb, iterations, maxNumCells)
if nargin < 3
    iterations = 150;
end
if nargin < 4
    maxNumCells = 10000;
end
if nargin < 2
    lb = 10;
end

ProgressBar(['Calculating gene-wise variation for dataset ''' ds.name ''''],true);


numCells = size(ds.data,2);
if (numCells > maxNumCells)
    ds = ds.randSample(maxNumCells);
    numCells = maxNumCells;%keep this line if we start using numCells below later
end

dstpm = TPM(mean(ds.data, 2));

sel = dstpm >= lb; 

ds = ds.geneSubset(sel);%filter genes to avoid problems with division by zero, etc, and to save compilation time. No point in looking at too lowly expressed ones anyway, they will not become significant.
numGenes = size(ds.data,1);


SNOLogCVS = zeros(numGenes, iterations);
SNOVariances = zeros(numGenes, iterations);
[logCVDS,varianceDS] = GetVarAndLogCV(ds);
%so we need to store the generated 
for it = 1:iterations
    SNO = GenSampDs(ds);
    [SNOLogCVS(:,it),SNOVariances(:,it)] = GetVarAndLogCV(SNO);
    ProgressBar(floor(100*it/iterations));
end

logCVSNOm = mean(SNOLogCVS,2);

logCVDifference = logCVDS - logCVSNOm;

genes = ds.genes;

%Calculate p value with a non-parametric method: Scale all genes to the
%same variance and translate so the mean is zero. Then use all genes for
%the non-parametric test, which in total gives a lot of data points.
%{
%scale and translate:
stddevvars = std(vvv, [], 2);
meanvars = mean(vvv, 2);
standardized = vvv - meanvars;
standardized = standardized ./ stddevvars;
%}


%Now calculate p values using a one sample z-test.
%This is an approximation, normality is rejected for 18% or so of the genes
%calculate the distribution from the SNO datasets

mn = mean(SNOLogCVS,2);
stddev = std(SNOLogCVS,0,2);

pVals = zeros(numGenes,1);
for g = 1:numGenes
    [~,pVals(g,1)] = ztest(logCVDS(g,1),mn(g,:),stddev(g,1),'Tail','right');
end

ProgressBar('Done');


end

function [logCV,variances] = GetVarAndLogCV(ds)
    ds = TPM(ds);
    avgRefExpr = mean(ds.data,2);
    variances = var(ds.data, 0, 2);
    sd = sqrt(variances);
    cv_ = sd ./ (avgRefExpr + 0.05);%Coefficient of Variation = std/mean. Adding 0.05, a neglectably small number, to handle too lowly expressed genes
    logCV = log2(cv_ + 1);%the + 1 says that no variance -> 0 value
end
 
%Generates a SNO dataset where gene expression is preserved instead of
%counts per cell
function ds = GenSampDs(d1)
    %create empty dataset
    ds = d1;
    numCells = size(d1.data,2);
    numGenes = size(ds.genes,1);
    ds.data = zeros(numGenes,numCells);%sparsify in the end


    %generate probabilities
    countsPerGene = sum(d1.data,2);
    countsPerCell = sum(d1.data,1);
    
    prob = countsPerCell./sum(countsPerCell,2);
    prob = full(prob);
    
    %generate a vector with the sum of probabilities up to this cell (including this cell)
    probSum = prob;%just create a vector of the right size
    probSum(1,1) = 0;%not really needed, but here for clarity
    for ii = 2:numCells
        probSum(1, ii) = probSum(1,ii-1) + prob(1, ii-1);
    end

    edges = [probSum max(probSum(1,end),1)];%add right edge; make sure it is not smaller than the previous due to roundoff problems. 
    edges = edges.';
        
    %generate genes
    for i = 1:numGenes
        if countsPerGene(i,1) > 0 %handle the fact that cells sometimes can contain 0 UMIs, which doesn't work with the code below. The data is automatically 0 in that case.
            %ProgressBar(floor(100*i/numCells));
            %generate n random numbers between 0-1 and sort them
            r = rand([countsPerGene(i,1) 1]);
            a = histcounts(r,edges);
            ds.data(i,:) = a.';
        end
    end

    ds.data = sparse(ds.data);

end


