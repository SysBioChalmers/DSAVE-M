function ret = DSAVEGetPrecalcTotVarFromBulk()
% DSAVEGetPrecalcTotVarFromBulk
%   Gets the DSAVE total variation for a few precalculated gene expression 
%   ranges. The data is generated from importable data if needed, but the
%   idea is that the user shall not be dependent on such data, and use the
%   .mat file only.
%
% Usage: bulkDataMeasurements = DSAVEGetPrecalcTotVarFromBulk();
%
% Johan Gustafsson, 2019-06-10
%

SCDep.init();
persistent v;
if isempty(v) 
    disp('reading DSAVE precalculated bulk total variation values...');
    prevDir = SCDep.setPathToSource();
    filename = '../data/DSAVE_bulk_tot_var.mat';
    if(~exist(filename,'file'))    
        disp('no .mat file found, generating data');
        
        %calculate mean deviance between bulk samples (CD4+ T cells) from Blueprint
        %read TMM normalized data
        samp = ImportTabSepSamples('../../ImportableData/tcellCD4ProfilesTMMNormalized.txt');

        bulkSamples = samp.sampleSubset(32:39);

        %Rescale the bulk samples so that they on average have 10^6 counts per
        %sample, even though the TMM normalization leads to that they will have a
        %different number of counts
        normFact = 10^6 / mean(sum(bulkSamples.data,1),2);
        bulkSamples.data = bulkSamples.data .* normFact;

        bulkMean1Vs1_05_100k = DSAVEGetTotalVariationFromBulk(bulkSamples, false, 100000, 0.5);
        bulkMean4Vs4_05_100k = DSAVEGetTotalVariationFromBulk(bulkSamples, true, 100000, 0.5);

        bulkMean1Vs1_05_2 = DSAVEGetTotalVariationFromBulk(bulkSamples, false, 2, 0.5);
        bulkMean4Vs4_05_2 = DSAVEGetTotalVariationFromBulk(bulkSamples, true, 2, 0.5);

        bulkMean1Vs1_100_100k = DSAVEGetTotalVariationFromBulk(bulkSamples, false, 100000, 100);
        bulkMean4Vs4_100_100k = DSAVEGetTotalVariationFromBulk(bulkSamples, true, 100000, 100);

        bulkMean1Vs1_05_50 = DSAVEGetTotalVariationFromBulk(bulkSamples, false, 50, 0.5);
        bulkMean4Vs4_05_50 = DSAVEGetTotalVariationFromBulk(bulkSamples, true, 50, 0.5);

        bulkMean1Vs1_2_100 = DSAVEGetTotalVariationFromBulk(bulkSamples, false, 100, 2);
        bulkMean4Vs4_2_100 = DSAVEGetTotalVariationFromBulk(bulkSamples, true, 100, 2);

        names = {'PseudoTPM: 0.5-100k', 'PseudoTPM: 0.5-2', 'PseudoTPM: 100-100k', ...
                 'PseudoTPM: 0.5-50', 'PseudoTPM: 2-100' };
        data_1vs1 = [bulkMean1Vs1_05_100k, bulkMean1Vs1_05_2, bulkMean1Vs1_100_100k, ...
                     bulkMean1Vs1_05_50, bulkMean1Vs1_2_100];
        data_4vs4 = [bulkMean4Vs4_05_100k, bulkMean4Vs4_05_2, bulkMean4Vs4_100_100k, ...
                     bulkMean4Vs4_05_50, bulkMean4Vs4_2_100];

        v.geneRanges = names;
        v.totalVariation_1vs1 = data_1vs1;
        v.totalVariation_4vs4 = data_4vs4;

        save(filename, 'v');
    else
        a = load(filename);
        v = a.v;
    end
    SCDep.restoreDir(prevDir);
end
ret = v;

end
