
iterations = 100000;
%iterations = 100;


%% Fig A

progbar = ProgrBar('Gene-wise: Fig A');


%breast cancer
bc2 = SCDep.scd_bc2;
bc2t = bc2.cellSubset(bc2.cellType == Celltype.TCellCD4Pos | bc2.cellType == Celltype.TCellCD8Pos | bc2.cellType == Celltype.TCellReg);
bc2tSub = bc2t.cellSubset(1:2000);%do not do random sampling for reproducability reasons.
[genesbc,logCVDifferencebc,pvalsbc,bcSNOVariances,bcSNOCountsPerGene] = DSAVEGetGeneVariation(bc2tSub,1,iterations,10000, progbar.GetSubContext(0.33));

%lung cancer
[lc,~] = SCDep.scd_lc;
lct = lc.cellSubset(lc.cellType == Celltype.TCellCD4Pos | lc.cellType == Celltype.TCellCD8Pos | lc.cellType == Celltype.TCellReg);
lctSub = lct.cellSubset(1:2000);%do not do random sampling for reproducability reasons.
[geneslc,logCVDifferencelc,pvalslc,lcSNOVariances,lcSNOCountsPerGene] = DSAVEGetGeneVariation(lctSub,1,iterations,10000, progbar.GetSubContext(0.33));

%hca
hcacb = SCDep.scd_hca_cb;
hca_cb1 = hcacb.cellSubset(strcmp(hcacb.sampleIds,'CB1'));
hcat = hca_cb1.cellSubset(hca_cb1.cellType == Celltype.TCell | hca_cb1.cellType == Celltype.TCellCD4Pos | hca_cb1.cellType == Celltype.TCellCD8Pos);
hcatSub = hcat.cellSubset(1:2000);%do not do random sampling for reproducability reasons.
[geneshca,logCVDifferencehca,pvalshca,hcaSNOVariances,hcaSNOCountsPerGene] = DSAVEGetGeneVariation(hcatSub,1,iterations,10000, progbar.GetSubContext(0.33));

progbar.Done();%finish progress bar here, some printouts below

pValsbcAdj = AdjustPval(pvalsbc,'benjamini',1);
sum(pValsbcAdj< 0.05)

[sortedValsBc,ibc] = sort(logCVDifferencebc.', 'descend');
xValsBc = 1:size(logCVDifferencebc,1);
sortedGenesBC = genesbc(ibc,:);
sortedPValsbcAdj = pValsbcAdj(ibc,:);

figure
plot(xValsBc(1,1:10:end),sortedValsBc(1,1:10:end),'LineWidth',2);%the plot looks bad if you have too many points, so we only use every 10th point
xlabel('Gene index')
ylabel('BTM variation')
title('Variation per Gene');
hold on

[sortedValsLc, ilc] = sort(logCVDifferencelc.', 'descend');
sortedGenesLC = geneslc(ilc,:);
xValsLc = 1:size(logCVDifferencelc,1);
plot(xValsLc(1,1:10:end),sortedValsLc(1,1:10:end),'LineWidth',2);
hold on
axis([0 14000 -0.1 2.5]);
pValslcAdj = AdjustPval(pvalslc,'benjamini',1);
sortedPValslcAdj = pValslcAdj(ilc,:);

[sortedValsHca, ihca] = sort(logCVDifferencehca.', 'descend');
sortedGenesHCA = geneshca(ihca,:);
xValsHca = 1:size(logCVDifferencehca,1);
plot(xValsHca(1,1:10:end),sortedValsHca(1,1:10:end),'LineWidth',2);
legend({'BC T cells','LC T cells', 'HCA CB T cells'});
set(gca,'FontSize',11);

pValshcaAdj = AdjustPval(pvalshca,'benjamini',1);
sortedPValshcaAdj = pValshcaAdj(ihca,:);



%some info: significant number of genes for the different datasets
sum(pValsbcAdj < 0.05)
sum(pValslcAdj < 0.05)
sum(pValshcaAdj < 0.05)

%% Fig 4A-D Supl. - Reproducability
progbar = ProgrBar('Gene-wise: Reproducibility');

%create a second run with 100,000 iterations
[genesbc2,logCVDifferencebc2,pvalsbc2,bcSNOVariances2,bcSNOCountsPerGene2] = DSAVEGetGeneVariation(bc2tSub,1,iterations,10000, progbar.GetSubContext(1));

progbar.Done();

%p value correlation between 2 runs - Fig B
figure
plot([0;10], [0;10])
hold on
scatter(pvalsbc, pvalsbc2);
xlabel('BTM Variation p-value run 1')
ylabel('BTM variation p-value run 2')
title('Gene-wise BTM Variation P-value Reproducibility, 100,000 It.');
axis([0 0.05 0 0.05]);
set(gca,'FontSize',11);

%for correlation - only look at the values below 0.05 - the rest is not
%interesting and will likely dominate the picture
pvalsbcbelow005 = pvalsbc(pvalsbc <= 0.05, :);
pvalsbc2below005 = pvalsbc2(pvalsbc <= 0.05, :);%select on the same dataset to get the same genes in both datasets
corrcoef(pvalsbcbelow005, pvalsbc2below005)

%run 2 runs with 100 iterations
[genesbc3,logCVDifferencebc3,pvalsbc3,bcSNOVariances3,bcSNOCountsPerGene3] = DSAVEGetGeneVariation(bc2tSub,1,100,10000);
[genesbc4,logCVDifferencebc4,pvalsbc4,bcSNOVariances4,bcSNOCountsPerGene4] = DSAVEGetGeneVariation(bc2tSub,1,100,10000);
%Fig A:
figure
plot([0;10], [0;10])
hold on
scatter(logCVDifferencebc3, logCVDifferencebc4);
xlabel('BTM Variation run 1')
ylabel('BTM variation run 2')
title('Gene-wise BTM Variation Reproducibility, 100 It.');
axis([0 3 0 3]);
set(gca,'FontSize',11);

corrcoef(logCVDifferencebc3, logCVDifferencebc4)


%run 2 runs with 10,000 iterations
[genesbc5,logCVDifferencebc5,pvalsbc5,bcSNOVariances5,bcSNOCountsPerGene5] = DSAVEGetGeneVariation(bc2tSub,1,10000,10000);
[genesbc6,logCVDifferencebc6,pvalsbc6,bcSNOVariances6,bcSNOCountsPerGene6] = DSAVEGetGeneVariation(bc2tSub,1,10000,10000);

%p value correlation between 2 runs
figure
plot([0;10], [0;10])
hold on
scatter(pvalsbc5, pvalsbc6);
xlabel('BTM Variation p-value run 1')
ylabel('BTM variation p-value run 2')
title('Gene-wise BTM Variation P-value Reproducibility, 10,000 It.');
axis([0 0.05 0 0.05]);
set(gca,'FontSize',11);

%for correlation - only look at the values below 0.05 - the rest is not
%interesting and will likely dominate the picture
pvalsbc5below005 = pvalsbc5(pvalsbc5 <= 0.05, :);
pvalsbc6below005 = pvalsbc6(pvalsbc5 <= 0.05, :);%select on the same dataset to get the same genes in both datasets
corrcoef(pvalsbc5below005, pvalsbc6below005)


%run 2 runs with 1,000 iterations
[genesbc7,logCVDifferencebc7,pvalsbc7,bcSNOVariances7,bcSNOCountsPerGene7] = DSAVEGetGeneVariation(bc2tSub,1,1000,10000);
[genesbc8,logCVDifferencebc8,pvalsbc8,bcSNOVariances8,bcSNOCountsPerGene8] = DSAVEGetGeneVariation(bc2tSub,1,1000,10000);

%p value correlation between 2 runs
figure
plot([0;10], [0;10])
hold on
scatter(pvalsbc7, pvalsbc8);
xlabel('BTM Variation p-value run 1')
ylabel('BTM variation p-value run 2')
title('Gene-wise BTM Variation P-value Reproducibility, 1,000 It.');
axis([0 0.05 0 0.05]);
set(gca,'FontSize',11);

%for correlation - only look at the values below 0.05 - the rest is not
%interesting and will likely dominate the picture
pvalsbc7below005 = pvalsbc7(pvalsbc7 <= 0.05, :);
pvalsbc8below005 = pvalsbc8(pvalsbc7 <= 0.05, :);%select on the same dataset to get the same genes in both datasets
corrcoef(pvalsbc7below005, pvalsbc8below005)



%% Fig D
%check what happens if we replace the 50/250 worst genes with SNO counterparts

templInfo = DSAVEGetStandardTemplate();

progbar = ProgrBar('Gene-wise: Fig D');

bc2tSub2 = bc2tSub;
lctSub2 = lctSub;
hcatSub2 = hcatSub;


%generate a template for each dataset
bcSNO = DSAVEGenerateSNODataset(bc2tSub2, progbar.GetSubContext(0.015));
lcSNO = DSAVEGenerateSNODataset(lctSub2, progbar.GetSubContext(0.015));
hcaSNO = DSAVEGenerateSNODataset(hcatSub2, progbar.GetSubContext(0.015));

templInfo2 = templInfo;
templInfo2.fractionUpperOutliers = 0;
templInfo2.fractionLowerOutliers = 0;
bcNoRemDSAVEScore = DSAVECalcBTMScore(bc2tSub2, templInfo2, progbar.GetSubContext(0.105));
lcNoRemDSAVEScore = DSAVECalcBTMScore(lctSub2, templInfo2, progbar.GetSubContext(0.105));
hcaNoRemDSAVEScore = DSAVECalcBTMScore(hcatSub2, templInfo2, progbar.GetSubContext(0.105));

bc2tSub2red = bc2tSub2;
worstGenesBc = genesbc(ibc(1:50));
[~,ia,~] = intersect(bc2tSub2red.genes,worstGenesBc);
bc2tSub2red.data(ia,:) = bcSNO.data(ia,:);
bcNoW50DSAVEScore = DSAVECalcBTMScore(bc2tSub2red, templInfo2, progbar.GetSubContext(0.105));

lctSub2red = lctSub2;
worstGenesLc = geneslc(ilc(1:50));
[~,ia,~] = intersect(lctSub2red.genes,worstGenesLc);
lctSub2red.data(ia,:) = lcSNO.data(ia,:);
lcNoW50DSAVEScore = DSAVECalcBTMScore(lctSub2red, templInfo2, progbar.GetSubContext(0.105));

hcatSub2red = hcatSub2;
worstGenesHca = geneshca(ihca(1:50));
[~,ia,~] = intersect(hcatSub2red.genes,worstGenesHca);
hcatSub2red.data(ia,:) = hcaSNO.data(ia,:);
hcaNoW50DSAVEScore = DSAVECalcBTMScore(hcatSub2red, templInfo2, progbar.GetSubContext(0.105));

bc2tSub2red2 = bc2tSub2;
worstGenesBc2 = genesbc(ibc(1:250));
[~,ia,~] = intersect(bc2tSub2red2.genes,worstGenesBc2);
bc2tSub2red2.data(ia,:) = bcSNO.data(ia,:);
bcNoW250DSAVEScore = DSAVECalcBTMScore(bc2tSub2red2, templInfo2, progbar.GetSubContext(0.105));

lctSub2red2 = lctSub2;
worstGenesLc2 = geneslc(ilc(1:250));
[~,ia,~] = intersect(lctSub2red2.genes,worstGenesLc2);
lctSub2red2.data(ia,:) = lcSNO.data(ia,:);
lcNoW250DSAVEScore = DSAVECalcBTMScore(lctSub2red2, templInfo2, progbar.GetSubContext(0.105));

hcatSub2red2 = hcatSub2;
worstGenesHca2 = geneshca(ihca(1:250));
[~,ia,~] = intersect(hcatSub2red2.genes,worstGenesHca2);
hcatSub2red2.data(ia,:) = hcaSNO.data(ia,:);
hcaNoW250DSAVEScore = DSAVECalcBTMScore(hcatSub2red2, templInfo2, progbar.GetSubContext(0.105));

progbar.Done();
disp('Copy the values below into the excel sheet, figure D.')
values = [bcNoRemDSAVEScore.DSAVEScore bcNoW50DSAVEScore.DSAVEScore bcNoW250DSAVEScore.DSAVEScore lcNoRemDSAVEScore.DSAVEScore lcNoW50DSAVEScore.DSAVEScore lcNoW250DSAVEScore.DSAVEScore hcaNoRemDSAVEScore.DSAVEScore hcaNoW50DSAVEScore.DSAVEScore hcaNoW250DSAVEScore.DSAVEScore]


%% Fig C - Wenn diagram over intersection of 50 most variable genes

worstGenesBCAdj = sortedGenesBC(sortedPValsbcAdj < 0.05);
worstGenesLCAdj = sortedGenesLC(sortedPValslcAdj < 0.05);
worstGenesHCAAdj = sortedGenesHCA(sortedPValshcaAdj < 0.05);
worstGenesBCAdj = worstGenesBCAdj(1:250,:);
worstGenesLCAdj = worstGenesLCAdj(1:250,:);
worstGenesHCAAdj = worstGenesHCAAdj(1:250,:);

[bc,lc,hca,bclc,bchca,lchca,all_] = CreateVennDiagramSets(worstGenesBCAdj,worstGenesLCAdj,worstGenesHCAAdj);

a = zeros(1,7);
a(1,1) = size(bc,1);
a(1,2) = size(lc,1);
a(1,3) = size(hca,1);
a(1,4) = size(bclc,1);
a(1,5) = size(bchca,1);
a(1,6) = size(lchca,1);
a(1,7) = size(all_,1);
a

bc
lc
hca
bclc
bchca
lchca
all_

%% Fig B Variation value at p = 0.05 (uncorrected) per gene expression

totUMIs = sum(sum(bc2tSub2.data,1),2);
meanExpr = 10^6*bcSNOCountsPerGene/totUMIs;

%sort the rows of the variances
varSorted = sort(bcSNOVariances,2);
lcvSorted = log(sqrt(varSorted) ./ (meanExpr+0.05) + 1);

numGenes = size(lcvSorted,1);
confvars005 = zeros(numGenes,1);
confvars001 = zeros(numGenes,1);
numRep = size(lcvSorted,2);

conf005Index = round(0.95 * numRep);
conf001Index = round(0.99 * numRep);

for g = 1:numGenes
    confvars005(g,1) = lcvSorted(g,conf005Index);
    %if there are several identical values in a row, this must be taken
    %care of if there is an identical value to the left; then we must find
    %the next value and take that
    if lcvSorted(g,conf005Index) == lcvSorted(g,conf005Index-1)
        nextIndex = conf005Index;
        while (nextIndex < numRep) && (lcvSorted(g,conf005Index) == lcvSorted(g,nextIndex))
            nextIndex = nextIndex + 1;
        end
        confvars005(g,1) = lcvSorted(g,nextIndex);
    end
    confvars001(g,1) = lcvSorted(g,conf001Index);
    if lcvSorted(g,conf001Index) == lcvSorted(g,conf001Index-1)
        nextIndex = conf001Index;
        while (nextIndex < numRep) && (lcvSorted(g,conf001Index) == lcvSorted(g,nextIndex))
            nextIndex = nextIndex + 1;
        end
        confvars001(g,1) = lcvSorted(g,nextIndex);    
    end
end

%calculate the BTM variation by subtracting the mean CV val (which is the
%sampling noise) from the value at the conf interval
meanCvs = mean(lcvSorted, 2);

btm005 = confvars005 - meanCvs;
btm001 = confvars001 - meanCvs;

figure
semilogx(meanExpr, btm001);
hold on
semilogx(meanExpr, btm005);
xlabel('Gene expression (CPM)')
ylabel('BTM variation')
title('Confidence Interval of BTM Variation per Gene Expression');
axis([0 5000 0 0.80]);
set(gca,'FontSize',11);
legend({'p = 0.01', 'p = 0.05'})

%Check that the results are reasonable:
cvsTest = log(sqrt(bcSNOVariances(4,:))/(meanExpr(4,:)+0.05) + 1);% TPM around 2
cvsMean = mean(cvsTest);
fakeBTM = cvsTest - cvsMean;
larg = fakeBTM >= btm005(4);
pval = sum(larg) / size(larg,2);%supposed to be 0.05, which it is
figure
histogram(cvsTest - cvsMean) 

