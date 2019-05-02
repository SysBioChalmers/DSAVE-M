function meanRes = EvaluateClusterSize(ds, maxSizeSet, upperBoundTPM, lowerBoundTPM)

ds_red = TPM(ds);

ProgressBar(strcat('Evaluate cluster size:',ds.name),true);

totset = ds_red.data;

%remove all genes not within expression range
avgRefExpr = mean(totset,2);
rem = ( avgRefExpr > upperBoundTPM ) | ( avgRefExpr < lowerBoundTPM );
totset(rem,:) = [];
totNumCells = size(totset,2);

steps = 100;
maxSize = min(maxSizeSet, totNumCells/2-1);
repetitions = 30;
res = zeros(2, steps*repetitions);
meanRes = zeros(2, steps);

%prepare the steps to produce a nice graph
numSampVector = zeros(100,1);
if (maxSize < 101)
    numSampVector(1:maxSize,1) = 1:maxSize;
    numSampVector(maxSize+1:100,1) = maxSize;
elseif (maxSize < 1000)
    step = floor((maxSize-1)/(steps-1));
    for i = 1:steps
        numSampVector(i) = 1 + (i-1)*step;
    end
else
    %make sure it is fine-grained enough below 500
    smallStep = 10;
    for i = 1:steps/2
        numSampVector(i) = 1 + (i-1)*smallStep;
    end
    start = numSampVector(steps/2);
    largeStep = floor((maxSize-start)/(steps/2));
    for i = steps/2+1:steps
        a = i-steps/2;
        numSampVector(i) = start + a*largeStep;
    end
end

for sz = 1:steps
    ProgressBar(sz/steps*100);
    numSamp = numSampVector(sz);%make sure the graphs start at 1.
    for repetition = 1:repetitions
        as = randsample(totNumCells, numSamp);
        a = mean(totset(:,as),2);
        ind = (1:totNumCells).';
        ind(as) = [];%remove the cells used in a to ensure we get other cells in b
        bsInd = randsample(totNumCells-numSamp, numSamp);
        bs = ind(bsInd);
        b = mean(totset(:,bs),2);
        
        geneScore = abs(log2((a+0.05)./(b+0.05)));
        index = (sz-1)*repetitions + repetition;
        res(1,index) = sz;
        res(2,index) = mean(geneScore);
    end
    meanRes(1, sz) = numSamp;
    st = ((sz-1)*repetitions + 1);
    en = ((sz-1)*repetitions+repetitions);
    meanRes(2, sz) = mean(res(2, st:en), 2);
end

ProgressBar('Done');

end
