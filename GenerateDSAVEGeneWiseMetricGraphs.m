
%% Fig A
%breast cancer
bc2 = SCDep.scd_bc2;
bc2t = bc2.cellSubset(bc2.paperClass == Celltype.TCellCD4Pos | bc2.paperClass == Celltype.TCellCD8Pos | bc2.paperClass == Celltype.TCellReg);
bc2tSub = bc2t.randSample(10000);
tic
[logCVDifferencebc,genesbc,bcSNOVariances, pvalsbc] = DSAVEGetGeneVariation(bc2tSub,1,150,10000);
t = toc

%figure
%histfit(bcSNOVariances(17,:));

%lung cancer
[lc,~] = SCDep.scd_lc;
lct = lc.cellSubset(lc.paperClass == Celltype.TCellCD4Pos | lc.paperClass == Celltype.TCellCD8Pos | lc.paperClass == Celltype.TCellReg);
lctSub = lct.randSample(10000);
tic
[logCVDifferencelc,geneslc,lcSNOVariances, pvalslc] = DSAVEGetGeneVariation(lctSub,1,150,10000);
t = toc

%hca
hcacb = SCDep.scd_hca_cb;
hca_cb1 = hcacb.cellSubset(strcmp(hcacb.sampleIds,'CB1'));
hcat = hca_cb1.cellSubset(hca_cb1.custClass == Celltype.TCell | hca_cb1.custClass == Celltype.TCellCD4Pos | hca_cb1.custClass == Celltype.TCellCD8Pos);
hcatSub = hcat.randSample(10000);
tic
[logCVDifferencehca,geneshca,hcaSNOVariances, pvalshca] = DSAVEGetGeneVariation(hcatSub,1,150,10000);
t = toc

%this code list genes that appear more than once
%[gg, ggg] = unique(hcacb.genes);
%aaa = hcacb.genes;
%aaa(ggg,:) = [];

%Save the data in case we want to recreate the graphs; it takes a whole day
%to run the lines above down to here...
prevDir = SCDep.setPathToSource();
save('../../TempData/logCVDifferencebc.mat','logCVDifferencebc');
save('../../TempData/genesbc.mat','genesbc');
save('../../TempData/pvalsbc.mat','pvalsbc');

save('../../TempData/logCVDifferencelc.mat','logCVDifferencelc');
save('../../TempData/geneslc.mat','geneslc');
save('../../TempData/pvalslc.mat','pvalslc');

save('../../TempData/logCVDifferencehca.mat','logCVDifferencehca');
save('../../TempData/geneshca.mat','geneshca');
save('../../TempData/pvalshca.mat','pvalshca');
SCDep.restoreDir(prevDir);

pValsbcAdj = AdjustPval(pvalsbc,'benjamini',1);
sum(pValsbcAdj< 0.05)
%sum(pvalsbc < 0.05)

[sortedValsBc,ibc] = sort(logCVDifferencebc.', 'descend');
xValsBc = 1:size(logCVDifferencebc,1);
figure
plot(xValsBc(1,1:10:end),sortedValsBc(1,1:10:end),'LineWidth',2);%the plot looks bad if you have too many points, so we only use every 10th point
xlabel('Gene index')
ylabel('BTM variation')
title('Variation per Gene');
hold on


[sortedValsLc, ilc] = sort(logCVDifferencelc.', 'descend');
xValsLc = 1:size(logCVDifferencelc,1);
plot(xValsLc(1,1:10:end),sortedValsLc(1,1:10:end),'LineWidth',2);
hold on
axis([0 14000 -0.1 2.5]);
pValslcAdj = AdjustPval(pvalslc,'benjamini',1);

[sortedValsHca, ihca] = sort(logCVDifferencehca.', 'descend');
xValsHca = 1:size(logCVDifferencehca,1);
plot(xValsHca(1,1:10:end),sortedValsHca(1,1:10:end),'LineWidth',2);
legend({'BC T cells','LC T cells', 'HCA CB T cells'});
set(gca,'FontSize',11);

pValshcaAdj = AdjustPval(pvalshca,'benjamini',1);



%some info: significant number of genes for the different datasets
sum(pValsbcAdj < 0.05)
sum(pValslcAdj < 0.05)
sum(pValshcaAdj < 0.05)

%

%% Corr Fig, not used
%{
%make correlation graph; doesn't look as good though
[~,ia,ib] = intersect(geneshca,geneslc);
figure
scatter(logCVDifferencehca(ia), logCVDifferencelc(ib));
xlabel('Bc T cells gene-wise BTM variation')
ylabel('Lc T cells gene-wise BTM variation')
title('Gene-wise BTM variation correlation between datasets.');
%}

%% Fig D

%check what happens if we replace the 50/250 worst genes with SNO counterparts
templInfo = DSAVEGetStandardTemplate();
bc2tSub2 = bc2tSub.randSample(2000);%make sure we compare the same datasets to avoid stochaisticity
lctSub2 = lctSub.randSample(2000);%make sure we compare the same datasets to avoid stochaisticity
hcatSub2 = hcatSub.randSample(2000);%make sure we compare the same datasets to avoid stochaisticity
bcDSAVEScore = CalcDSAVE(bc2tSub2, templInfo);
lcDSAVEScore = CalcDSAVE(lctSub2, templInfo);
hcaDSAVEScore = CalcDSAVE(hcatSub2, templInfo);

%generate a template for each dataset
bcSNO = GenerateSamplingSSDataset(bc2tSub2);
lcSNO = GenerateSamplingSSDataset(lctSub2);
hcaSNO = GenerateSamplingSSDataset(hcatSub2);

templInfo2 = templInfo;
templInfo2.fractionUpperOutliers = 0;
templInfo2.fractionLowerOutliers = 0;
bcNoRemDSAVEScore = CalcDSAVE(bc2tSub2, templInfo2);
lcNoRemDSAVEScore = CalcDSAVE(lctSub2, templInfo2);
hcaNoRemDSAVEScore = CalcDSAVE(hcatSub2, templInfo2);

bc2tSub2red = bc2tSub2;
worstGenesBc = genesbc(ibc(1:50));
[~,ia,~] = intersect(bc2tSub2red.genes,worstGenesBc);
bc2tSub2red.data(ia,:) = bcSNO.data(ia,:);
bcNoW50DSAVEScore = CalcDSAVE(bc2tSub2red, templInfo2);

lctSub2red = lctSub2;
worstGenesLc = geneslc(ilc(1:50));
[~,ia,~] = intersect(lctSub2red.genes,worstGenesLc);
lctSub2red.data(ia,:) = lcSNO.data(ia,:);
lcNoW50DSAVEScore = CalcDSAVE(lctSub2red, templInfo2);

hcatSub2red = hcatSub2;
worstGenesHca = geneshca(ihca(1:50));
[~,ia,~] = intersect(hcatSub2red.genes,worstGenesHca);
hcatSub2red.data(ia,:) = hcaSNO.data(ia,:);
hcaNoW50DSAVEScore = CalcDSAVE(hcatSub2red, templInfo2);

bc2tSub2red2 = bc2tSub2;
worstGenesBc2 = genesbc(ibc(1:250));
[~,ia,~] = intersect(bc2tSub2red2.genes,worstGenesBc2);
bc2tSub2red2.data(ia,:) = bcSNO.data(ia,:);
bcNoW250DSAVEScore = CalcDSAVE(bc2tSub2red2, templInfo2);

lctSub2red2 = lctSub2;
worstGenesLc2 = geneslc(ilc(1:250));
[~,ia,~] = intersect(lctSub2red2.genes,worstGenesLc2);
lctSub2red2.data(ia,:) = lcSNO.data(ia,:);
lcNoW250DSAVEScore = CalcDSAVE(lctSub2red2, templInfo2);

hcatSub2red2 = hcatSub2;
worstGenesHca2 = geneshca(ihca(1:250));
[~,ia,~] = intersect(hcatSub2red2.genes,worstGenesHca2);
hcatSub2red2.data(ia,:) = hcaSNO.data(ia,:);
hcaNoW250DSAVEScore = CalcDSAVE(hcatSub2red2, templInfo2);

values = [bcNoRemDSAVEScore.DSAVEScore bcNoW50DSAVEScore.DSAVEScore bcNoW250DSAVEScore.DSAVEScore lcNoRemDSAVEScore.DSAVEScore lcNoW50DSAVEScore.DSAVEScore lcNoW250DSAVEScore.DSAVEScore hcaNoRemDSAVEScore.DSAVEScore hcaNoW50DSAVEScore.DSAVEScore hcaNoW250DSAVEScore.DSAVEScore]

%% Test if the most variable genes are mainly lowly expressed:
%[logCVDifferencebc,genesbc,bcSNOVariances, pvalsbc] = DSAVEGetGeneVariation(bc2tSub,1,150,10000);
%[sortedValsBc,ibc] = sort(logCVDifferencebc.', 'descend');

[~,ia,~] = intersect(bc2tSub.genes,worstGenesBc);
meanexprbc = mean(TPM(bc2tSub.data),2);
meanexprbclowly = full(meanexprbc(ia,:));
data50worst = full(bc2tSub.data(ia,:));
a = sum(data50worst > 1,2);
format longG
b = [meanexprbclowly a]
format short % default
%conclusion: Most of the most variably genes are lowly expressed; but the
%variation comes from a few cells only (i.e. cells with more than one count). Therefore filter on 10 CPM, since those are not that interesting.


%% Fig C - Wenn diagram over intersection of 50 most variable genes

%get the worst genes above 10 CPM
%worstGenesBc = genesbc(ibc(1:50));
%sortedGenesBc = genesbc(ibc);
%{
%first synchronize the genes between the dataset and the results from the
%variation calculation:
meanExprBc = mean(TPM(bc2tSub.data),2);
[commonGenesBc,ia,ib] = intersect(bc2tSub.genes,genesbc);
commonMeanExprBc = meanExprBc(ia,:);
commonLogCVDifferenceBc = logCVDifferencebc(ib,:);
%Then filter out all genes below 10 CPM:
above10 = commonMeanExprBc >= 10;
commonGenesAbove10Bc = commonGenesBc(above10,:);
commonMeanExprAbove10Bc = commonMeanExprBc(above10,:);
commonLogCVDifferenceAbove10Bc = commonLogCVDifferenceBc(above10,:);
%now sort according to variation
[sortedValsAbove10Bc,iAbove10Bc] = sort(commonLogCVDifferenceAbove10Bc.', 'descend');
commonGenesAbove10Bc = commonGenesAbove10Bc(iAbove10Bc);
commonGenesAbove10_50WorstBc = commonGenesAbove10Bc(1:250,:);

meanExprLc = mean(TPM(lctSub.data),2);
[commonGenesLc,ia,ib] = intersect(lctSub.genes,geneslc);
commonMeanExprLc = meanExprLc(ia,:);
commonLogCVDifferenceLc = logCVDifferencelc(ib,:);
%Then filter out all genes below 10 CPM:
above10 = commonMeanExprLc >= 10;
commonGenesAbove10Lc = commonGenesLc(above10,:);
commonMeanExprAbove10Lc = commonMeanExprLc(above10,:);
commonLogCVDifferenceAbove10Lc = commonLogCVDifferenceLc(above10,:);
%now sort according to variation
[sortedValsAbove10Lc,iAbove10Lc] = sort(commonLogCVDifferenceAbove10Lc.', 'descend');
commonGenesAbove10Lc = commonGenesAbove10Lc(iAbove10Lc);
commonGenesAbove10_50WorstLc = commonGenesAbove10Lc(1:250,:);

meanExprHca = mean(TPM(hcatSub.data),2);
[commonGenesHca,ia,ib] = intersect(hcatSub.genes,geneshca);
commonMeanExprHca = meanExprHca(ia,:);
commonLogCVDifferenceHca = logCVDifferencehca(ib,:);
%Then filter out all genes below 10 CPM:
above10 = commonMeanExprHca >= 10;
commonGenesAbove10Hca = commonGenesHca(above10,:);
commonMeanExprAbove10Hca = commonMeanExprHca(above10,:);
commonLogCVDifferenceAbove10Hca = commonLogCVDifferenceHca(above10,:);
%now sort according to variation
[sortedValsAbove10Hca,iAbove10Hca] = sort(commonLogCVDifferenceAbove10Hca.', 'descend');
commonGenesAbove10Hca = commonGenesAbove10Hca(iAbove10Hca);
commonGenesAbove10_50WorstHca = commonGenesAbove10Hca(1:250,:);



%logCVDifferencebc(strcmp(genesbc,'MALAT1'))
%meanExprBc(strcmp(bc2tSub.genes,'MALAT1'))
logCVDifferencebc(strcmp(genesbc,'HBB'))
meanExprBc(strcmp(bc2tSub.genes,'HBB'))
dataHBB = bc2tSub.data(strcmp(bc2tSub.genes,'HBB'),:);
sel = dataHBB > 0;
dataHBB(sel)

dataHBB = lctSub.data(strcmp(lctSub.genes,'HBB'),:);
sel = dataHBB > 0;
dataHBB(sel)
figure
histogram(dataHBB(sel))

dataHBB = hcatSub.data(strcmp(hcatSub.genes,'HBB'),:);
sel = dataHBB > 0;
dataHBB(sel)
figure
histogram(dataHBB(sel))

dataGNLY = lctSub.data(strcmp(lctSub.genes,'GNLY'),:);
sel = dataGNLY > 0;
dataGNLY(sel)
figure
histogram(dataGNLY(sel))


%numSorted = size(sortedGenesBc,1);
%sortIndex = 1:numSorted;

%}

%[bc,lc,hca,bclc,bchca,lchca,all] = CreateVennDiagramSets(worstGenesBc,worstGenesLc,worstGenesHca);
[bc,lc,hca,bclc,bchca,lchca,all] = CreateVennDiagramSets(worstGenesBc2,worstGenesLc2,worstGenesHca2);
%[bc,lc,hca,bclc,bchca,lchca,all] = CreateVennDiagramSets(commonGenesAbove10_50WorstBc,commonGenesAbove10_50WorstLc,commonGenesAbove10_50WorstHca);
%{
bc_lc = intersect(worstGenesBc,worstGenesLc);
bc_hca = intersect(worstGenesBc,worstGenesHca);
lc_hca = intersect(worstGenesLc,worstGenesHca);

bc_lc_hca = intersect(bc_lc,worstGenesHca);
%}
a = zeros(1,7);
a(1,1) = size(bc,1);
a(1,2) = size(lc,1);
a(1,3) = size(hca,1);
a(1,4) = size(bclc,1);
a(1,5) = size(bchca,1);
a(1,6) = size(lchca,1);
a(1,7) = size(all,1);
a

bc
lc
hca
bclc
bchca
lchca
all

%% Fig B Variation value at p = 0.05 (uncorrected) per gene expression

%we don't have all the genes in the results
[~,ia,ib] = intersect(bc2tSub2.genes,genesbc);

%first calculate mean expression
meanExpr = mean(TPM(bc2tSub2.data(ia,:)),2);

variances = var(bcSNOVariances,0,2);

%need to fix the indices of the variances as well to match the intersect
%above
variances = variances(ib,:);

numGenes = size(genesbc,1);
confvars = zeros(numGenes,1);
means = zeros(numGenes,1);

for g = 1:numGenes
    %fit variances to a normal distribution
    pd = fitdist(bcSNOVariances(g,:).','Normal');
    ci = paramci(pd,'Alpha',.1);%for one-sided test
    confvars(g,1) = ci(2,1);%this is the value at the confidence interval; the value means the variance of the data at the conf interval
    means(g,1) = pd.mu;
end

%recalculate to log2(cv + 1);
confcvs = (confvars + 0.05) ./ (meanExpr+0.05);
logconfcvs = log2(confcvs + 1);

meancvs = means ./ meanExpr;
logmeancvs = log2(meancvs + 1);

diffs = logconfcvs - logmeancvs;

%figure
%scatter(meanExpr,diffs);

%skip for the genes < 1 TPM
sel = meanExpr >= 1 & meanExpr <= 5000;

dfs = diffs(sel);
exp = meanExpr(sel);
%now sort the values
[es,ie] = sort(exp);
rs = dfs(ie);

%create a plot using sliding window of 300 genes (makes it somewhat smooth)
numPoints = size(rs,1) - 299;
xs = zeros(1,numPoints);
ys = zeros(1,numPoints);

for i = 1:numPoints
    ys(1,i) = mean(rs(i:i+299,1));
    xs(1,i) = mean(es(i:i+299,1));
end

figure
h = area(xs,ys,'LineWidth',1);
xlabel('Gene expression (CPM)')
ylabel('BTM variation')
title('Confidence Interval of BTM Variation per Gene Expression');
axis([0 100 0 0.04]);
h(1).FaceColor = [0.75 0.75 1]
set(gca,'FontSize',11);


figure
histfit(log2(bcSNOVariances(94,:)));


%% The Anderson-Darling test for normality on the variances
sel = meanExpr >= 1;

varsel = log2(sqrt(bcSNOVariances(sel,:)));
numGenes = size(varsel,1);
normPVals = zeros(numGenes, 1);
for i = 1:numGenes
    [~,normPVals(i,1)] = adtest(varsel(i,:));
end
sum(normPVals < 0.05)
perc = 1-sum(normPVals < 0.05)/numGenes
rejPerc = sum(normPVals < 0.05)/numGenes
