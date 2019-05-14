%Calculates the average pairwise total variation between a list of bulk samples, with a TPM
%filtration on the genes. If pool4samples is true, it will compare the mean
%of 4 randomly selected samples with the mean of 4 others. This expects a
%list of 8 samples.
function res = DSAVEGetTotalVariationFromBulk(s, pool4samples, upperBoundTPM, lowerBoundTPM)
    means = mean(s.data,2);
    badGenes = means < lowerBoundTPM | means > upperBoundTPM;
    s = s.geneSubset(s.genes(~badGenes)); 
    numSamp = size(s.sampleIds,2);
    numGenes = size(s.genes,1);
    diffs = zeros(numGenes,numSamp*(numSamp-1)/2);
    index = 1;

    if ~pool4samples
        for i = 1:numSamp-1
           for j = i+1:numSamp
               diffs(:,index) = log((s.data(:,i)+0.05)./(s.data(:,j)+0.05));
               index = index + 1;
           end
        end
    else
        if size(s.data,2) ~= 8
            error('4 on 4 only works with 8 samples');
        end
        %here the average of 4 samples is compared to the
        %average of another 4.
        %get all permutations of four
        ind = 1:numSamp;
        k = round(numSamp/2);
        combs = nchoosek(ind,k);
        numCombs = size(combs,1);
        diffs = zeros(numGenes,numCombs);
        for i = 1:numCombs
            a = s.data(:,combs(i,:));
            notind = 1:numSamp;
            notind(combs(i,:)) = [];
            b = s.data(:,notind);
            ma = mean(a,2);
            mb = mean(b,2);
            diffs(:,i) = log((ma+0.05)./(mb+0.05));
        end
    end
    res = mean(mean(abs(diffs),2),1);
end
