
%% Fig A

progbar = ProgrBar('Gene-wise: Fig A');


%breast cancer
bc2 = SCDep.scd_bc2;
bc2t = bc2.cellSubset(bc2.cellType == Celltype.TCellCD4Pos | bc2.cellType == Celltype.TCellCD8Pos | bc2.cellType == Celltype.TCellReg);
bc2tSub = bc2t.randSample(10000);
[logCVDifferencebc,genesbc,bcSNOVariances, pvalsbc] = DSAVEGetGeneVariation(bc2tSub,1,150,10000, progbar.GetSubContext(0.33));

%figure
%histfit(bcSNOVariances(17,:));

%lung cancer
[lc,~] = SCDep.scd_lc;
lct = lc.cellSubset(lc.cellType == Celltype.TCellCD4Pos | lc.cellType == Celltype.TCellCD8Pos | lc.cellType == Celltype.TCellReg);
lctSub = lct.randSample(10000);
[logCVDifferencelc,geneslc,lcSNOVariances, pvalslc] = DSAVEGetGeneVariation(lctSub,1,150,10000, progbar.GetSubContext(0.33));

%hca
hcacb = SCDep.scd_hca_cb;
hca_cb1 = hcacb.cellSubset(strcmp(hcacb.sampleIds,'CB1'));
hcat = hca_cb1.cellSubset(hca_cb1.cellType == Celltype.TCell | hca_cb1.cellType == Celltype.TCellCD4Pos | hca_cb1.cellType == Celltype.TCellCD8Pos);
hcatSub = hcat.randSample(10000);
[logCVDifferencehca,geneshca,hcaSNOVariances, pvalshca] = DSAVEGetGeneVariation(hcatSub,1,150,10000, progbar.GetSubContext(0.33));

progbar.Done();%finish progress bar here, some printouts below

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

progbar = ProgrBar('Gene-wise: Fig D');

bc2tSub2 = bc2tSub.randSample(2000);%make sure we compare the same datasets to avoid stochaisticity
lctSub2 = lctSub.randSample(2000);%make sure we compare the same datasets to avoid stochaisticity
hcatSub2 = hcatSub.randSample(2000);%make sure we compare the same datasets to avoid stochaisticity

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

[bc,lc,hca,bclc,bchca,lchca,all_] = CreateVennDiagramSets(worstGenesBc2,worstGenesLc2,worstGenesHca2);

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
exp_ = meanExpr(sel);
%now sort the values
[es,ie] = sort(exp_);
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

%% some tests

bcsubs = bc2tSub.geneSubset(genesbc);%this works since geneSubset is deterministic; however very shaky, the order of the genes could differ
mns = mean(TPM(bcsubs.data),2);
counts = sum(bcsubs.data,2);
vvv = log2(sqrt(bcSNOVariances)./(mns+0.05) + 1);%this is the same logCV as calculated for the observed
www = log(sqrt(bcSNOVariances)./(mns+0.00000005));


%try to overlay them, scaled by variance
stddevvars = std(vvv, [], 2);
meanvars = mean(vvv, 2);
standardized = vvv - meanvars;
standardized = standardized ./ stddevvars;
stddevvars2 = std(www, [], 2);
meanvars2 = mean(www, 2);
standardized2 = www - meanvars2;
standardized2 = standardized2 ./ stddevvars2;
figure
histfit(standardized(82,:));

edges = -6:0.5:6;
xes = -5.75:0.5:5.75


X = [ones(length(counts),1) log(counts)];

%hockeyklubban:

%startläge:
figure
scatter(log(counts),www(:,1))

cc = www+.5*log(counts);%.5 samma som sqrt på counts - detta transformerar CV så att lutningen mot counts nästan försvinner
figure
scatter(log(counts),cc(:,1))

aa = cc-mean(cc,2);%centrera, gör att cv-värdena ligger runt x-axeln (och fixar till den lilla uträtning som inte skedde i förra steget)
bb = www - mean(www,2);%bb is identical to aa, showing that the first step can be skipped if the second step is performed

figure
scatter(log(counts),aa(:,1))

figure
scatter(log(counts),bb(:,1))


X = [ones(length(counts),1) log(counts)];
b = X\(aa.^2)

stddev = std(aa,[],2);
figure
scatter(log(counts), aa(:,1)./stddev);%div med stddev gör att spridning blir någorlunda jämn

figure
scatter(log(counts), log(stddev))

figure
scatter(log(mns), log(stddev))

rnge2 = 1:0.1:6;
rnge = 10.^rnge2;
prob = rnge ./ 10^6;
numUMI = 1500;
%mean = n*p, var = n*p(1-p)
cvggg = sqrt(1-prob) ./ sqrt(numUMI.*prob);
figure
scatter(log(rnge), log(cvggg))


%for i = 1:150
%    scatter(log(counts), aa(:,i)./stddev);%div med stddev gör att spridning blir jämn, men introducerar en lutning igen
%end

%detta är själva grafen där man ser att höguttryckta gener avviker från
%linjen
stddev = std(aa,[],2);
figure
scatter(-0.5*log(counts), log(stddev))



X = [ones(length(counts),1) log(counts)];
b = X\(aa.^2)


%%%%%%%% HERE
aa = www+0.5*log(counts);
% correct mean for log(count) effect
% BUT fudge factor plays a role (TPM) same for all genes so subtract
% overall mean
aa = aa-mean(mean(aa));
figure
for i = 1:150 
 scatter(log(counts),aa(:,i))
 hold on
end
%%% std(log(CV)) is directly proportional to -.5*log(counts)
stddev = std(aa,[],2);
figure
for i = 1:150 
 scatter(log(counts),aa(:,i)./stddev)
 hold on
end


%generate a dataset with the same number of UMIs in each cell
hcasame = SCDep.scd_hca_cb.cellSubset(1:30000);
UMI1 = sum(SCDep.scd_hca_cb.data,1);
hcasame = SCDep.scd_hca_cb.cellSubset(UMI1 > 3000);
templInfoSpec = DSAVEGetStandardTemplate();
templInfoSpec.UMIDistr = repmat(2000,1,2000);
hcasame = DSAVEAlignDataset(hcasame, templInfoSpec);
templInfoSpec2 = templInfoSpec;
templInfoSpec2.UMIDistr = repmat(1000,1,2000);
hcasame2 = DSAVEAlignDataset(hcasame, templInfoSpec2);


%SNO = DSAVEGenerateSNODataset(hcasame);

[logCVDifferenceSame,genesSame,sameSNOVariances, pvalssame] = DSAVEGetGeneVariation(hcasame,1,150,10000);
[logCVDifferenceSame2,genesSame2,sameSNOVariances2, pvalssame2] = DSAVEGetGeneVariation(hcasame2,1,150,10000);

SNOsub = hcasame.geneSubset(genesSame);%this works since geneSubset is deterministic; however very shaky, the order of the genes could differ
mnssame = mean(TPM(SNOsub.data),2);
countssame = sum(SNOsub.data,2);
vvvsame = log2(sqrt(sameSNOVariances)./(mnssame+0.05) + 1);%this is the same logCV as calculated for the observed
wwwsame = log(sqrt(sameSNOVariances)./(mnssame+0.00000005));

SNOsub2 = hcasame2.geneSubset(genesSame2);%this works since geneSubset is deterministic; however very shaky, the order of the genes could differ
mnssame2 = mean(TPM(SNOsub2.data),2);
countssame2 = sum(SNOsub2.data,2);
vvvsame2 = log2(sqrt(sameSNOVariances2)./(mnssame2+0.05) + 1);%this is the same logCV as calculated for the observed
wwwsame2 = log(sqrt(sameSNOVariances2)./(mnssame2+0.00000005));

%test that gene synch worked
all(strcmp(SNOsub.genes, genesSame)) %yes it did!


ccsame = wwwsame+.5*log(countssame);%.5 samma som sqrt på counts - detta transformerar CV så att lutningen mot counts nästan försvinner
aasame = ccsame-mean(ccsame,2);%centrera, gör att cv-värdena ligger runt x-axeln (och fixar till den lilla uträtning som inte skedde i förra steget)

figure
scatter(log(countssame),ccsame(:,3))


figure
scatter(log(countssame),aasame(:,3))

figure
scatter(log(counts),bb(:,1))


X = [ones(length(counts),1) log(counts)];
b = X\(aa.^2)

stddevsame = std(aasame,[],2);
figure
scatter(log(countssame), aasame(:,1)./stddevsame);%div med stddev gör att spridning blir någorlunda jämn

figure
scatter(log(countssame), log(stddevsame))

figure
scatter(log(counts), log(stddev))

rnge2 = 1:0.1:6;
rnge = 10.^rnge2;
prob = rnge ./ 10^6;
numUMI = 1500;
%mean = n*p, var = n*p(1-p)
cvggg = sqrt(1-prob) ./ sqrt(numUMI.*prob);
figure
scatter(log(rnge), log(cvggg))

%detta är själva grafen där man ser att höguttryckta gener avviker från
%linjen
stddev = std(aa,[],2);
figure
scatter(-0.5*log(counts), log(stddev))






%kör med buckets
stddevvarssame = std(wwwsame, [], 2);
meanvarssame = mean(wwwsame, 2);
standardizedsame = wwwsame - meanvarssame;
standardizedsame = standardizedsame ./ stddevvarssame;

stddevvarssame2 = std(wwwsame2, [], 2);
meanvarssame2 = mean(wwwsame2, 2);
standardizedsame2 = wwwsame2 - meanvarssame2;
standardizedsame2 = standardizedsame2 ./ stddevvarssame2;

%loop through a whole bunch of buckets
stnddata = standardizedsame;
lowerBounds = 1:20;
upperBounds = 2:21;
lowerBounds = 100:500:1100;
upperBounds = [600 1100 10000000];

lowerBounds = 50:10:100;
upperBounds = lowerBounds + 10;

lowerBounds = 10:20;
upperBounds = lowerBounds + 1;

%plot histogram:
figure
for i = 1:size(lowerBounds,2)
%for i = 1:3
    sel = mnssame <= upperBounds(1,i) & mnssame >= lowerBounds(1,i);
    dat = stnddata(sel,:);
    [aa,ab] = size(dat);
    aB = reshape(dat, [1, aa*ab]);
    edges = -12:0.25:12;
    %xes = -3.75:0.5:3.75
    xes = -11.875:0.25:11.875

    ac = histcounts(aB,edges);
    ac = ac / size(aB,2);
    plot(xes,ac)
    hold on
   
end

%plot integrated histogram:
figure
%for i = 1:size(lowerBounds,2)
for i = 1:size(lowerBounds,2)
    sel = mnssame <= upperBounds(1,i) & mnssame >= lowerBounds(1,i);
    dat = stnddata(sel,:);
    [aa,ab] = size(dat);
    aB = reshape(dat, [1, aa*ab]);
    edges = -12:0.25:12;
    %xes = -3.75:0.5:3.75
    xes = -11.875:0.25:11.875;

    ac = histcounts(aB,edges);
    ac = ac / size(aB,2);
    accum = cumsum(ac);
    plot(xes,accum)
    hold on
   
end
legend(arrayfun(@num2str,lowerBounds,'UniformOutput',false));

%plot integrated histogram for 1000 UMI set:
figure
%for i = 1:size(lowerBounds,2)
for i = 1:size(lowerBounds,2)
    sel = mnssame2 <= upperBounds(1,i) & mnssame2 >= lowerBounds(1,i);
    dat = stnddata(sel,:);
    [aa,ab] = size(dat);
    aB = reshape(dat, [1, aa*ab]);
    edges = -12:0.25:12;
    %xes = -3.75:0.5:3.75
    xes = -11.875:0.25:11.875;

    ac = histcounts(aB,edges);
    ac = ac / size(aB,2);
    accum = cumsum(ac);
    plot(xes,accum)
    hold on
   
end
legend(arrayfun(@num2str,lowerBounds,'UniformOutput',false));



%put the genes into buckets and look at the standardized distributions for each bucket:
stnddata = standardizedsame;
sel05_2 = mnssame < 2 & mnssame > 0.5;
dat05_2a = stnddata(sel05_2,:);
[aa,ab] = size(dat05_2a);
aB05_2 = reshape(dat05_2a, [1, aa*ab]);
edges = -6:0.25:6;
%xes = -3.75:0.5:3.75
xes = -5.875:0.25:5.875

figure
ac = histcounts(aB05_2,edges);
ac = ac / size(aB05_2,2);
plot(xes,ac)
hold on

sel5_10 = mnssame < 10 & mnssame > 5;
dat5_10a = stnddata(sel5_10,:);
[aa,ab] = size(dat5_10a);
aB5_10 = reshape(dat5_10a, [1, aa*ab]);
ac = histcounts(aB5_10,edges);
ac = ac / size(aB5_10,2);
plot(xes,ac)
hold on

sel20_100 = mnssame < 100 & mnssame > 20;
dat20_100a = stnddata(sel20_100,:);
[aa,ab] = size(dat20_100a);
aB20_100 = reshape(dat20_100a, [1, aa*ab]);
ac = histcounts(aB20_100,edges);
ac = ac / size(aB20_100,2);
plot(xes,ac)
hold on

selL100 = mnssame > 100;
datL100a = stnddata(selL100,:);
[aa,ab] = size(datL100a);
aBL100 = reshape(datL100a, [1, aa*ab]);
ac = histcounts(aBL100,edges);
ac = ac / size(aBL100,2);
plot(xes,ac)
legend({'0.5-2 TPM','5-10 TPM','20-100 TPM','>100 TPM'});










hold on

sort05_2 = sort(aB05_2);
ind = round(size(aB05_2,2)*0.05);
plot([aB05_2(1,aB05_2),aB05_2(1,aB05_2)], [0,0.1]);










figure
for i = 1:size(standardized,1)
    a = histcounts(standardized(i,:),edges);
    plot(xes,a)
    hold on
end

figure
histfit(stddevvars)

%rescale and translate observed values with std and mean from SNO
%first reverse engineer the total variation:
%hmm, this is actually not needed if translation is done with mean; since translation will exactly undo this!
SNOMean = mean(vvv,2);
totalVar = logCVDifferencebc + SNOMean;
%then translate and rescale (where translation is not needed)
obsStandardized = logCVDifferencebc ./ stddevvars;
%plot the real values
figure
histfit(obsStandardized.');

%Now calc p values
%First pool all values for all genes into one vector
[a,b] = size(standardized)
B = reshape(standardized, [1, a*b]);
newPVals = zeros(b,1);
for i = 1:b
    numLower = sum(B < obsStandardized(i,1));
    numHigher = sum(B > obsStandardized(i,1));
    newPVals(i,1) = numHigher / (numLower + numHigher);%ignore all values that are the same
end
[newPVals(1:100,1) pvalsbc(1:100,1)]

%takes a really long time
%c = zeros(a,1);
%for i = 1:a
%    [~,c(i,1)] = kstest2(B, standardized2(i,:));
%end
%{
d = c(1:351,1);
sum(d < 0.05) / size(d,1)
figure
histogram(d)

h = zeros(a,1);
C = B(:,randsample(size(B,2),15000));

size(unique(C))
figure
histogram(C,150)
%}
%{
for i = 1:a
    i
    [~,h(i,1)] = kstest2(C, standardized2(i,:));
end

figure
histogram(h,50)
%}

e = zeros(1000,1);
f = randperm(1000);
for i = 1:1000
    disp(i);
    e(i,1) = ranksum(standardized(i+1,:), standardized(i,:));

    %[~,e(i,1)] = kstest2(standardized2(i+1,:), standardized2(i,:));
    %[~,e(i,1)] = kstest2(chi2rnd(5, 1, 150), chi2rnd(5, 1, 150));
    
end

figure
histogram(e)

%plot CVs against total gene count
figure
for i = 1:5
    scatter(counts, vvv(:,i));
    hold on
end

figure
for i = 1:20
    scatter(counts, 2.^ vvv(:,i) - 1);
    hold on
end


figure
for i = 1:20
    scatter(log2(counts), vvv(:,i));
    hold on
end




figure
scatter(log2(sqrt(counts)), sqrt(vvv(:,i)));

figure
scatter(log2(sqrt(counts)), www(:,1));


figure
histogram(standardized2(17,:), 50);
size(unique(standardized2(17,:)))
figure
histogram(e, 50);


g = zeros(1000,1);
for i = 1:1000
    disp(i);
    %[~,e(i,1)] = kstest2(standardized2(f(1,i),:), standardized2(i,:));
    [~,g(i,1)] = kstest2(chi2rnd(5, 1, 10000), chi2rnd(5, 1, 150));
    
end
figure
histogram(g,100)


%put the genes into buckets and look at the standardized distributions for each bucket:
stnddata = standardized2;
sel05_2 = mns < 2 & mns > 0.5;
dat05_2a = stnddata(sel05_2,:);
[aa,ab] = size(dat05_2a);
aB05_2 = reshape(dat05_2a, [1, aa*ab]);
edges = -6:0.25:6;
%xes = -3.75:0.5:3.75
xes = -5.875:0.25:5.875

figure
ac = histcounts(aB05_2,edges);
ac = ac / size(aB05_2,2);
plot(xes,ac)
hold on

sel5_10 = mns < 10 & mns > 5;
dat5_10a = stnddata(sel5_10,:);
[aa,ab] = size(dat5_10a);
aB5_10 = reshape(dat5_10a, [1, aa*ab]);
ac = histcounts(aB5_10,edges);
ac = ac / size(aB5_10,2);
plot(xes,ac)
hold on

sel20_100 = mns < 100 & mns > 20;
dat20_100a = stnddata(sel20_100,:);
[aa,ab] = size(dat20_100a);
aB20_100 = reshape(dat20_100a, [1, aa*ab]);
ac = histcounts(aB20_100,edges);
ac = ac / size(aB20_100,2);
plot(xes,ac)
hold on

selL100 = mns > 100;
datL100a = stnddata(selL100,:);
[aa,ab] = size(datL100a);
aBL100 = reshape(datL100a, [1, aa*ab]);
ac = histcounts(aBL100,edges);
ac = ac / size(aBL100,2);
plot(xes,ac)
legend({'0.5-2 TPM','5-10 TPM','20-100 TPM','>100 TPM'});

hold on
sort05_2 = sort(aB05_2);
ind = round(size(aB05_2,2)*0.05);
plot([aB05_2(1,aB05_2),aB05_2(1,aB05_2)], [0,0.1]);

figure
histfit(vvv(4,:).');
hold on
plot([totalVar(4,1), totalVar(4,1)], [0,20]);

figure
histfit(standardized2(4,:).');
hold on
plot([obsStandardized(4,1), obsStandardized(4,1)], [0,20])

figure
histfit(vvv(5,:).');
hold on
plot([totalVar(5,1), totalVar(5,1)], [0,20]);

figure
histfit(standardized2(5,:).');
hold on
plot([obsStandardized(5,1), obsStandardized(5,1)], [0,20])

figure
histfit(vvv(6,:).');
hold on
plot([totalVar(6,1), totalVar(6,1)], [0,20]);

figure
histfit(standardized2(6,:).');
hold on
plot([obsStandardized(6,1), obsStandardized(6,1)], [0,20])

%Now try to match a Xsq (Chi Square) curve, by first making sure min of
%each gene is 0 (i.e. translate), then scale by dividing with the mean:
chiish = vvv - min(vvv,[],2);
chiish = chiish ./ mean(chiish,2);
figure
histfit(chiish(85,:));
hold on
xesxsq = 0:0.1:5;
b = chi2pdf(xesxsq, 11)*150/5;
plot(xesxsq * 0.18, b*150/5)

%test to fit a gamma distribution:
edges2 = 0:0.25:4;
edges2xes = 0.125:0.25:3.875;
c = histcounts(chiish(85,:),edges2);
[phat,pci] = mle(chiish(85,:),'distribution','gamma');
d = gampdf(xesxsq,phat(:,1), phat(:,2));
figure
plot(edges2xes, c);
hold on
plot(xesxsq,d*150*0.25);


figure
for i = 1:100
    a = histcounts(chiish(i,:),edges2);
    plot(xes,a)
    hold on
end


%histcounts(X,nbins)


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
