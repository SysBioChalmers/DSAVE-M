function results = CalcDSAVE(origDs, templInfo, skipAlignment)

if nargin < 3
    skipAlignment = false;
end

if skipAlignment
   disp('Warning: skipping alignment!'); 
end

iterations = 15;
%cellPoolSize = 5000;

lbnonlog = 10;
ubnonlog = 1000;
numPoints = 1000;%use the same number of points to simplify for data allocation
meansLog = linspace(log10(lbnonlog), log10(ubnonlog), numPoints);
xes = 10 .^ meansLog;


numXVals = size(xes,2);

alignedCVs = zeros(iterations,numXVals);
samplingCVs = zeros(iterations,numXVals);
differenceCVs = zeros(iterations,numXVals);

for it = 1:iterations
    disp(strcat('running iteration: ', num2str(it)));
    if (skipAlignment)
        aligned = origDs;
    else
        aligned = DSAVEAlignDataset(origDs, templInfo);%use a different subset of cells in each loop iteration
    end
    alignedGeneCVs = GetGeneCVs(aligned, templInfo);
    
    SNO = GenerateSamplingSSDataset(aligned, size(aligned.data, 2));
    SNOGeneCVs = GetGeneCVs(SNO, templInfo);
    
    %throw away the most and least variable genes
    numGenes = size(aligned.data,1);
    numToDiscardUpper = round(templInfo.fractionUpperOutliers * numGenes);
    numToDiscardLower = round(templInfo.fractionLowerOutliers * numGenes);
    difference = alignedGeneCVs - SNOGeneCVs;
    [~,gi] = sort(difference);
    discard = [gi(1:numToDiscardLower); gi((end-numToDiscardUpper+1):end)];
    anyToDiscard = sum(discard,1) > 0;
    
    alignedGeneCVsRem = alignedGeneCVs;
    SNOGeneCVsRem = SNOGeneCVs;
    if anyToDiscard
        alignedGeneCVsRem(discard,:) = [];
        SNOGeneCVsRem(discard,:) = [];
    end
    
    alData = TPM(aligned.data);
    alData(discard,:) = [];
    SNOData = TPM(SNO.data);
    SNOData(discard,:) = [];
    
    
    [alCVs, alXes] = AverageIntoBins(alData, alignedGeneCVsRem, templInfo);
    [saCVs, saXes] = AverageIntoBins(SNOData, SNOGeneCVsRem, templInfo);

    %remove any NaNs (let linear interpolation fill in)
    
    %now, recalculate the x:es and cvs to the Xes in the template using linear interpolation;
    %otherwise it will be difficult to take the median later
    
    %first remove any points with duplicate x; these will
    %otherwise mess up the linear interpolation.
    [alXes,ia,~] = unique(alXes);
    alCVs = alCVs(1,ia);
    [saXes,ia,~] = unique(saXes);
    saCVs = saCVs(1,ia);
    %then use linear interpolation
    %if the code fails on either of these two lines, that is because there are no genes that fits some bins
    %in that case, try using a template that discards fewer outliers
    alignedCVs(it,:) = interp1(alXes,alCVs,xes);
    samplingCVs(it,:) = interp1(saXes,saCVs,xes);
    
    differenceCVs(it,:) = alignedCVs(it,:) - samplingCVs(it,:);%add a neglectable number to avoid any potential division by 0.
end

%now, take the median of the iterations from all three curves. 

results = struct();
results.alignedCVs = median(alignedCVs,1);
results.samplingCVs = median(samplingCVs,1);
results.differenceCVs = median(differenceCVs,1);
results.DSAVEScore = mean(results.differenceCVs);%mean over all points ranging over different TPM
results.tpms = xes;

end


function logcv = GetGeneCVs(ds, templInfo)
ds_red = TPM(ds);

%ProgressBar(strcat('Calculating logCV for dataset ''', ds.name, ''' as function of gene expr'),true);

totset = ds_red.data;
numGenes = size(totset, 1);

avgRefExpr = mean(totset,2);
%{
%precreate the random selection of subsets and calculate the means of those
numRepetitions = templInfo.cellPoolIterations;
poolData = zeros(numGenes, numRepetitions);
for repetition = 1:numRepetitions
    ProgressBar(repetition/numRepetitions*100);
    as = datasample(totset, templInfo.cellPoolSize, 2);%get numSamp number of cells with replacement (meaning one cell can be taken several times)
    poolData(:,repetition) = TPM(sum(as,2));%take the mean over the cells -> a vertical vector with the mean value for each gene
end
%}
%calc variances
variances = var(totset, 0, 2);
sd = sqrt(variances);
cv_ = sd ./ (avgRefExpr + 0.05);%Coefficient of Variation = std/mean. Adding 0.05, a neglectably small number, to handle too lowly expressed genes
logcv = log2(cv_ + 1);%the + 1 says that no variance -> 0 value

%ProgressBar('Done');

end



function [cv, meanGeneExpr] = AverageIntoBins(TPMData, logcv, templInfo)

avgRefExpr = mean(TPMData,2);


numBins = size(templInfo.binningInfo.means,2);

cv = zeros(1,numBins);
meanGeneExpr = templInfo.binningInfo.means;


for i = 1:numBins
    %select the genes within the expression range
    sel = avgRefExpr >= templInfo.binningInfo.lbs(1,i) & avgRefExpr <= templInfo.binningInfo.ubs(1,i);
    cv(1, i) = mean(logcv(sel));%y value in the graph
    meanGeneExpr(1,i) = 2^mean(log2(avgRefExpr(sel)+0.05)) - 0.05;%geometric-ish mean
end


end


