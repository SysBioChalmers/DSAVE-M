    
%% Fig A and B
%breast cancer
bc2 = SCDep.scd_bc2;
bc2t = bc2.cellSubset(bc2.paperClass == Celltype.TCellCD4Pos | bc2.paperClass == Celltype.TCellCD8Pos | bc2.paperClass == Celltype.TCellReg);
bc2t_bc2ln = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC2_LYMPHNODE'));
bc2tSub = bc2t_bc2ln.randSample(2500);
llsbc2 = DSAVEGetSingleCellDivergence(bc2tSub, 200);

%plot against number of umis
%numUMIs = sum(bc2tSub.data,1);
[bcx,bci] = sort(llsbc2);
%bcy = numUMIs(1,bci);
%figure
%scatter(bcx,bcy);

templInfo = DSAVEGetStandardTemplate();

%get DSAVE score for the original dataset and for the one where the 500
%most divergent have been removed
bcDSAVEScore = CalcDSAVE(bc2tSub, templInfo);

%now remove the 500 worst cells
bc2tSub500Less = bc2tSub.cellSubset(bci(1,501:2500));
bcDSAVEScore500Less = CalcDSAVE(bc2tSub500Less, templInfo);
%{
lls500Less = DSAVEGetSingleCellDivergence(bc2tSub500Less, 200);
numUMIs500Less = sum(bc2tSub500Less.data,1);
[bcx2,bci2] = sort(lls500Less);
bcy2 = numUMIs500Less(1,bci2);
figure
scatter(bcx2,bcy2);
%}
%now for the SNO
bcSNO = GenerateSamplingSSDataset(bc2tSub);
llsbcSNO = DSAVEGetSingleCellDivergence(bcSNO, 200);
[bcsnox,bcsnoi] = sort(llsbcSNO);


%lc
[lc,~] = SCDep.scd_lc;
lct = lc.cellSubset(lc.paperClass == Celltype.TCellCD4Pos | lc.paperClass == Celltype.TCellCD8Pos | lc.paperClass == Celltype.TCellReg);
lctSub2 = lct.randSample(2500);
llslc2 = DSAVEGetSingleCellDivergence(lctSub2, 200);
[lcx,lci] = sort(llslc2);
%sno
lcSNO = GenerateSamplingSSDataset(lctSub2);
llslcSNO = DSAVEGetSingleCellDivergence(lcSNO, 200);
[lcsnox,lcsnoi] = sort(llslcSNO);

%get DSAVE score for the original dataset and for the one where the 500
%most divergent have been removed
lctDSAVEScore = CalcDSAVE(lctSub2, templInfo);

%now remove the 500 worst cells
lct500Less = lctSub2.cellSubset(lci(1,501:2500));
lctDSAVEScore500Less = CalcDSAVE(lct500Less, templInfo);



%hca
hcacb = SCDep.scd_hca_cb;
hca_cb1 = hcacb.cellSubset(strcmp(hcacb.sampleIds,'CB1'));
hcat = hca_cb1.cellSubset(hca_cb1.custClass == Celltype.TCell | hca_cb1.custClass == Celltype.TCellCD4Pos | hca_cb1.custClass == Celltype.TCellCD8Pos);
hcab = hca_cb1.cellSubset(hca_cb1.custClass == Celltype.BCell);
hcam = hca_cb1.cellSubset(hca_cb1.custClass == Celltype.Monocyte);
hcatSub2 = hcat.randSample(2500);


llshcat = DSAVEGetSingleCellDivergence(hcatSub2, 200);
hcatDSAVEScore = CalcDSAVE(hcatSub2, templInfo);
[hcatx,hcati] = sort(llshcat);
hcat500Less = hcatSub2.cellSubset(hcati(1,501:2500));
llshcat500Less = DSAVEGetSingleCellDivergence(hcat500Less, 200);
%{
numUMIshcat500Less = sum(hcat500Less.data,1);
figure
%scatter(hcamixx,hcamixy,3,colorIndices);
scatter(llshcat500Less,numUMIshcat500Less,3);
%}
hcat500LessDSAVEScore = CalcDSAVE(hcat500Less, templInfo);

hcaSNO = GenerateSamplingSSDataset(hcatSub2);
llshcaSNO = DSAVEGetSingleCellDivergence(hcaSNO, 200);
[hcatsnox,hcatsnoi] = sort(llshcaSNO);


%B10k
b10k = SCDep.scd_pbmcb10000;
b10kSub = b10k.randSample(2500);
llsb10k = DSAVEGetSingleCellDivergence(b10kSub, 200);
[b10kx,b10ki] = sort(llsb10k);
%sno
b10kSNO = GenerateSamplingSSDataset(b10kSub);
llsb10kSNO = DSAVEGetSingleCellDivergence(b10kSNO, 200);
[b10ksnox,b10ksnoi] = sort(llsb10kSNO);

%get DSAVE score for the original dataset and for the one where the 500
%most divergent have been removed
b10kDSAVEScore = CalcDSAVE(b10kSub, templInfo);

%now remove the 500 worst cells
b10k500Less = b10kSub.cellSubset(b10ki(1,501:2500));
b10kDSAVEScore500Less = CalcDSAVE(b10k500Less, templInfo);


%Now plot fig A:
xs = 1:2500;
figure
plot(xs, -bcx, '-', 'Color', [0, 0.4470, 0.7410],'LineWidth',2);
hold on
plot(xs, -bcsnox, '--', 'Color', [0, 0.4470, 0.7410],'LineWidth',2);
hold on
%plot(xs, -lcx, '-', 'Color', [0.8500, 0.3250, 0.0980],'LineWidth',2);
plot(xs, -b10kx, '-', 'Color', [0.8500, 0.3250, 0.0980],'LineWidth',2);
hold on
%plot(xs, -lcsnox, '--', 'Color', [0.8500, 0.3250, 0.0980],'LineWidth',2);
plot(xs, -b10ksnox, '--', 'Color', [0.8500, 0.3250, 0.0980],'LineWidth',2);
hold on
plot(xs, -hcatx, '-', 'Color',[0.800, 0.600, 0.100],'LineWidth',2);
hold on
plot(xs, -hcatsnox, '--', 'Color',[0.800, 0.600, 0.100],'LineWidth',2);%[0.9290, 0.6940, 0.1250]
hold on
xlabel('Cell index')
ylabel('-Log likelihood')
title('Cell Divergence');
legend({'BC LN T cells, single pat', 'BC LN T cells, single pat - SNO', 'B10k B cells, single pat', 'B10k B cells, single pat - SNO', 'HCA CB T cells, single pat', 'HCA CB T cells, single pat - SNO'});
%axis([0 2500 500 2500]);
set(gca,'FontSize',11);

%Data for Fig B:
disp('Fig 5B data, copy to excel:');
a = [bcDSAVEScore.DSAVEScore bcDSAVEScore500Less.DSAVEScore b10kDSAVEScore.DSAVEScore b10kDSAVEScore500Less.DSAVEScore hcatDSAVEScore.DSAVEScore hcat500LessDSAVEScore.DSAVEScore]

%% Investigate why some datasets get a higher BTM variation after removing the most divergent cells
%{
%Tried BC4_BLOOD, gave a reduction, ok: 1.0772 -> 1.0648
%It seems bc4_tumor is more problematic:
bc2t_bc4blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
bc2tbc4bSub = bc2t_bc4blood.randSample(3500);
llsbc2bc4b = DSAVEGetSingleCellDivergence(bc2tbc4bSub, 200);

bc2t_bc4blood5 = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_BLOOD'));
bc2tbc4bSub5 = bc2t_bc4blood5.randSample(3500);
llsbc2bc4b5 = DSAVEGetSingleCellDivergence(bc2tbc4bSub5, 200);
[bcxbc4b5,bcibc4b5] = sort(llsbc2bc4b5);


%plot against number of umis
%numUMIs = sum(bc2tSub.data,1);
[bcxbc4b,bcibc4b] = sort(llsbc2bc4b);
%bcy = numUMIs(1,bci);
%figure
%scatter(bcx,bcy);

templInfo = DSAVEGetStandardTemplate();


%loop through and calculate the DSAVE score by removing x number of worst
%cells up to 1500
maxRem = 1500;
pts = 20
rvals = zeros(1,pts+1);
torems = zeros(1,pts+1);
for i = 0:pts
    toRem = i/pts * maxRem;
    subset = bc2tbc4bSub.cellSubset(bcibc4b(1,(toRem+1):3500));
    sc = CalcDSAVE(subset, templInfo);
    rvals(1,i+1) = sc.DSAVEScore;
    torems(1,i+1) = toRem;
end
figure
plot(torems,rvals);

val0 = CalcDSAVE(bc2tbc4bSub, templInfo);
bc2tbc4bSub1500less = bc2tbc4bSub.cellSubset(bcibc4b(1,1501:3500));
val1500 = CalcDSAVE(bc2tbc4bSub1500less, templInfo);
figure
plot(val0.tpms,val0.differenceCVs);
hold on
plot(val1500.tpms,val1500.differenceCVs);
mean0 = TPM(mean(bc2tbc4bSub.data,2));
geneslarg10_0 = sum(mean0 >= 10);
mean1500 = TPM(mean(bc2tbc4bSub1500less.data,2));
geneslarg10_1500 = sum(mean1500 >= 10)

%calc DSAVE on the genes that are lost:
sel = (mean0 >= 10) & (mean1500 < 10);
bc2tbc4bSubLostGenes = bc2tbc4bSub.geneSubset(sel);
[logCVDifference,gs,poolData1,poolDataS, pvals] = DSAVEGetGeneVariation(bc2tbc4bSubLostGenes,1,500,150,10000);
bc2tbc4bSubRestOfGenes = bc2tbc4bSub.geneSubset((~sel) & (mean0 >= 10));
[logCVDifference2,gs2,poolData12,poolDataS2, pvals2] = DSAVEGetGeneVariation(bc2tbc4bSubRestOfGenes,1,500,150,10000);
mean(logCVDifference)
meanslost = mean0(sel);
meanslost1500 = mean1500(sel);
[mean(logCVDifference) mean(logCVDifference2)]
%need to filter the logCVDifference2 genes since they were not filtered for >10
%TPM
%genesLarger10 = bc2tbc4bSub.genes(mean0 >= 10);
%[~,ia,~] = intersect(bc2tbc4bSubRestOfGenes.genes, genesLarger10);

%mean(logCVDifference2(ia))

%for blood
maxRem = 1500;
pts = 20;
rvals5 = zeros(1,pts+1);
torems5 = zeros(1,pts+1);
for i = 0:pts
    toRem = i/pts * maxRem;
    subset = bc2tbc4bSub5.cellSubset(bcibc4b5(1,(toRem+1):3500));
    sc = CalcDSAVE(subset, templInfo);
    rvals5(1,i+1) = sc.DSAVEScore;
    torems5(1,i+1) = toRem;
end
figure
plot(torems5,rvals5);

pcaData = full(bc2tbc4bSub.data(:,bcibc4b).');%unsparse and sort
%change order of data to make the red be painted over the greenish
pcaData = [pcaData(1501:3500,:);pcaData(1:1500,:)];
[coeff, score, latent] = pca(pcaData);
colors = [repmat(1,2000,1);repmat(0,1500,1)];
figure
gscatter(score(:,1),score(:,2),colors);

pcaData2 = full(bc2tbc4bSub5.data(:,bcibc4b5).');%unsparse and sort
%change order of data to make the red be painted over the greenish
pcaData2 = [pcaData2(1501:3500,:);pcaData2(1:1500,:)];
[coeff2, score2, latent2] = pca(pcaData2);
colors = [repmat(1,2000,1);repmat(0,1500,1)];
figure
gscatter(score2(:,1),score2(:,2),colors);


%get DSAVE score for the original dataset and for the one where the 500
%most divergent have been removed
bc2tbc4bDSAVEScore = CalcDSAVE(bc2tbc4bSub, templInfo);

%now remove the 500 worst cells
bc2tbc4bSub500Less = bc2tbc4bSub.cellSubset(bcibc4b(1,501:2500));
bc2tbc4bDSAVEScore500Less = CalcDSAVE(bc2tbc4bSub500Less, templInfo);
[bc2tbc4bDSAVEScore.DSAVEScore bc2tbc4bDSAVEScore500Less.DSAVEScore]

llsbc2bc4b500Less = DSAVEGetSingleCellDivergence(bc2tbc4bSub500Less, 200);
[bcxbc4b500Less,bcibc4b500Less] = sort(llsbc2bc4b500Less);
figure
plot(bcxbc4b)
hold on
plot([501:2500],bcxbc4b500Less)
figure
plot(bc2tbc4bDSAVEScore.alignedCVs)
hold on
plot(bc2tbc4bDSAVEScore500Less.alignedCVs)

figure
plot(bc2tbc4bDSAVEScore.samplingCVs)
hold on
plot(bc2tbc4bDSAVEScore500Less.samplingCVs)

figure
plot(bc2tbc4bDSAVEScore.differenceCVs)
hold on
plot(bc2tbc4bDSAVEScore500Less.differenceCVs)

%compare to the other bc2 patient
figure
plot(bcDSAVEScore.differenceCVs)
hold on
plot(bcDSAVEScore500Less.differenceCVs)

%experiment with different lower CPM bound
llsbc2bc4b2 = DSAVEGetSingleCellDivergence(bc2tbc4bSub, 200,100);

%plot against number of umis
%numUMIs = sum(bc2tSub.data,1);
[bcxbc4b2,bcibc4b2] = sort(llsbc2bc4b2);
%bcy = numUMIs(1,bci);
%figure
%scatter(bcx,bcy);

templInfo = DSAVEGetStandardTemplate();

%now remove the 500 worst cells
bc2tbc4bSub500Less2 = bc2tbc4bSub.cellSubset(bcibc4b2(1,501:2500));
bc2tbc4bDSAVEScore500Less2 = CalcDSAVE(bc2tbc4bSub500Less2, templInfo);
[bc2tbc4bDSAVEScore.DSAVEScore bc2tbc4bDSAVEScore500Less2.DSAVEScore]

umicounts = sum(bc2tbc4bSub.data,1);
mean(umicounts)
mean(umicounts(bcibc4b(1,501:2500)))

pcaData = full(bc2tbc4bSub.data(:,bcibc4b).');%unsparse and sort
%change order of data to make the red be painted over the greenish
pcaData = [pcaData(1001:2500,:);pcaData(1:1000,:)];
[coeff, score, latent] = pca(pcaData);
colors = [repmat(1,1500,1);repmat(0,1000,1)];
figure
gscatter(score(:,1),score(:,2),colors);




pcaDatahca = full(hcatSub2.data(:,hcati).');%unsparse and sort
[coeffhca, scorehca, latenthca] = pca(pcaDatahca);
colors = [repmat(0,500,1);repmat(1,2000,1)];
figure
gscatter(scorehca(:,1),scorehca(:,2),colors);

uc = sum(hcatSub2.data,1);
figure
histogram(uc)
llshcat2 = DSAVEGetSingleCellDivergence(hcatSub2, 1000);
[hcatx2,hcati2] = sort(llshcat2);
hcat500Less2 = hcatSub2.cellSubset(hcati2(1,501:2500));
hcat500LessDSAVEScore2 = CalcDSAVE(hcat500Less2, templInfo);
[hcatDSAVEScore.DSAVEScore hcat500LessDSAVEScore.DSAVEScore]

pcaDatahca2 = full(hcatSub2.data(:,hcati2).');%unsparse and sort
%change order of data to make the red be painted over the greenish
pcaDatahca2 = [pcaDatahca2(501:2500,:);pcaDatahca2(1:500,:)];
[coeffhca2, scorehca2, latenthca2] = pca(pcaDatahca2);
colors = [repmat(1,2000,1);repmat(0,500,1)];
figure
gscatter(scorehca2(:,1),scorehca2(:,2),colors);

%experiment with different mixes of two cell types to see if a 50% mix
%gives the maximum variation:
hcat2500 = hcat.randSample(2500);
hcam2500 = hcam.randSample(2500);
mix0 = hcat2500;
mix5 = hcat2500.randSample(2375).innerJoin(hcam2500.randSample(125));
mix10 = hcat2500.randSample(2250).innerJoin(hcam2500.randSample(250));
mix15 = hcat2500.randSample(2125).innerJoin(hcam2500.randSample(375));
mix20 = hcat2500.randSample(2000).innerJoin(hcam2500.randSample(500));
mix25 = hcat2500.randSample(1875).innerJoin(hcam2500.randSample(625));
mix30 = hcat2500.randSample(1750).innerJoin(hcam2500.randSample(750));
mix35 = hcat2500.randSample(1625).innerJoin(hcam2500.randSample(875));
mix40 = hcat2500.randSample(1500).innerJoin(hcam2500.randSample(1000));
mix45 = hcat2500.randSample(1375).innerJoin(hcam2500.randSample(1125));
mix50 = hcat2500.randSample(1250).innerJoin(hcam2500.randSample(1250));
score0 = CalcDSAVE(mix0, templInfo);
score5 = CalcDSAVE(mix5, templInfo);
score10 = CalcDSAVE(mix10, templInfo);
score15 = CalcDSAVE(mix15, templInfo);
score20 = CalcDSAVE(mix20, templInfo);
score25 = CalcDSAVE(mix25, templInfo);
score30 = CalcDSAVE(mix30, templInfo);
score35 = CalcDSAVE(mix35, templInfo);
score40 = CalcDSAVE(mix40, templInfo);
score45 = CalcDSAVE(mix45, templInfo);
score50 = CalcDSAVE(mix50, templInfo);
vals = [score0.DSAVEScore score5.DSAVEScore score10.DSAVEScore score15.DSAVEScore score20.DSAVEScore score25.DSAVEScore score30.DSAVEScore score35.DSAVEScore score40.DSAVEScore score45.DSAVEScore score50.DSAVEScore]
xes = [0 5 10 15 20 25 30 35 40 45 50];
figure
plot(xes,vals);
xlabel('Percent monocytes in mix');
ylabel('DSAVE variation score');
%}

%% Fig C - Num UMIs vs log likelihood
numUMIshcatSub2 = sum(hcatSub2.data,1);

linevalYs = 1500:500:10500;
linevalXes = zeros(1,size(linevalYs,2));
%generate line
for i = 1:size(linevalYs,2);
   ii = linevalYs(i);
   lb = ii - 250;
   ub = ii + 250;
   sel = numUMIshcatSub2 >= lb & numUMIshcatSub2 <= ub;
   linevalXes(1,i) = mean(llshcat(sel));
end



figure
scatter(llshcat,numUMIshcatSub2,3);
hold on
plot(linevalXes,linevalYs,'Color',[0, 0.8, 0.0],'LineWidth',2);
xlabel('Log likelihood')
ylabel('UMI counts')
title('UMI Counts vs Cell Divergence');
legend({'Individual cell', 'Mean log-likelihood'});
set(gca,'FontSize',11);
axis([-4500 -1400 0 14000]);

%% Fig D. MT-genes

matches = regexp(bc2tSub.genes, '^(MT-.*)$', 'match');
a = cellfun(@isempty, matches);
b = ~a;
sum(b)
bc2tSub.genes(b)

sumMito = sum(bc2tSub.data(b,:),1);
sumUMIs = sum(bc2tSub.data,1);
percMito = sumMito ./ sumUMIs;
%{
%Sliding window: this doesn't look as good
[sortpm, ia] = sort(percMito);
sortlls = llsbc2(1,ia);
numsum = 50;
numpoints = size(llsbc2,2) - numsum + 1;
linevalXes = zeros(1,numpoints);
linevalYs = zeros(1,numpoints);

for i = 1:numpoints
    linevalXes(1,i) = mean(sortlls(1,i:i+numsum-1),2);
    linevalYs(1,i) = mean(sortpm(1,i:i+numsum-1),2);
end
%}
linevalYs = 0.01:0.01:0.2;
linevalXes = zeros(1,20);
%generate line
for i = 1:20;
   ii = i*.01;
   lb = ii - 0.01;
   ub = ii + 0.01;
   sel = percMito >= lb & percMito <= ub;
   linevalXes(1,i) = mean(llsbc2(sel));
end


figure
scatter(llsbc2,percMito,3);
hold on
plot(linevalXes,linevalYs,'Color',[0, 0.8, 0.0],'LineWidth',2);

xlabel('Log likelihood')
ylabel('Fraction mitochondrial gene counts')
title('Fraction MT Counts vs Cell Divergence');
legend({'Individual cell', 'Mean log-likelihood'});
set(gca,'FontSize',11);

%% Fig 4B in supplementary - Number of detected genes
detGenes = sum(hcatSub2.data > 0,1);

linevalYs = 500:50:2400;
linevalXes = zeros(1,size(linevalYs,2));
%generate line
for i = 1:size(linevalYs,2);
   ii = linevalYs(i);
   lb = ii - 250;
   ub = ii + 250;
   sel = detGenes >= lb & detGenes <= ub;
   linevalXes(1,i) = mean(llshcat(sel));
end



figure
scatter(llshcat,detGenes,3);
hold on
plot(linevalXes,linevalYs,'Color',[0, 0.8, 0.0],'LineWidth',2);
xlabel('Log likelihood')
ylabel('Number of detected genes')
title('Number of Detected Genes vs Cell Divergence');
legend({'Individual cell', 'Mean log-likelihood'});
set(gca,'FontSize',11);
axis([-4500 -1400 450 2800]);


%% Fig F - Mix of B and T cells

%{
%This is code for generating the old plot where we showed a mix of b and t
cells
hcat2500 = hcat.randSample(2500);
hcam2500 = hcam.randSample(2500);

%hcabSub = hcab.randSample(500);
hcatSub = hcat2500.randSample(2250);
hcamSub = hcam2500.randSample(250);
hcamix = hcatSub.innerJoin(hcamSub);%the 500 first cells will be b cells
llshcamix = DSAVEGetSingleCellDivergence(hcamix, 200);

colorIndices = [repmat([1;0;0],1,2250) repmat([0;0;1],1,250)];

%plot against number of umis
numUMIs = sum(hcamix.data,1);


figure
scatter(llshcamix(1,1:2250),numUMIs(1,1:2250),3,colorIndices(:,1:2250).');
hold on
scatter(llshcamix(1,2251:2500),numUMIs(1,2251:2500),3,colorIndices(:,2251:2500).');
xlabel('Log likelihood')
ylabel('UMI counts')
title('UMI Counts vs Cell Divergence for a Mix of Cell Types');
legend({'T cells', 'Monocytes'});
axis([-2400 -1000 0 14000]);
set(gca,'FontSize',11);

hcamSub2 = hcam2500.randSample(2250);
hcatSub2 = hcat2500.randSample(250);
hcamix2 = hcamSub2.innerJoin(hcatSub2);%the 250 last cells will be t cells
llshcamix2 = DSAVEGetSingleCellDivergence(hcamix2, 200);

colorIndices = [repmat([1;0;0],1,2250) repmat([0;0;1],1,250)];

%plot against number of umis
numUMIs = sum(hcamix2.data,1);


figure
scatter(llshcamix2(1,1:2250),numUMIs(1,1:2250),3,colorIndices(:,1:2250).');
hold on
scatter(llshcamix2(1,2251:2500),numUMIs(1,2251:2500),3,colorIndices(:,2251:2500).');
xlabel('Log likelihood')
ylabel('UMI counts')
title('UMI Counts vs Cell Divergence for a Mix of B and T Cells');
legend({'T cells', 'B cells'});
axis([-2400 -1000 0 14000]);
set(gca,'FontSize',11);

%get pure t cells and monocytes as well
llshcapureT = DSAVEGetSingleCellDivergence(hcat2500, 200);
llshcapureM = DSAVEGetSingleCellDivergence(hcam2500, 200);

%test boxplot
%generate the data for the boxplot
boxdata = [llshcapureT llshcapureM llshcamix llshcamix2(:,2251:2500) llshcamix2(:,1:2250)].';
groupdata = [repmat(0,2500,1);repmat(1,2500,1);repmat(2,2250,1);repmat(3,250,1);repmat(4,2250,1);repmat(5,250,1)];

figure
boxplot(boxdata, groupdata, 'Labels', {'100% T cells' '100% Mon.' '90% T Cells','10% Mon.', '10% T Cells','90% Mon.'});
set(gca,'FontSize',11);
%}





%% Fig 4A in supplementary. Technical verification - run1 vs run 2.
%run divergence 2 times and plot

hcat2500 = hcat.randSample(2500);
hcam2500 = hcam.randSample(2500);

hcatSub = hcat2500.randSample(2250);
hcamSub = hcam2500.randSample(250);
hcamix = hcatSub.innerJoin(hcamSub);%the 500 first cells will be b cells


llshcamix1 = DSAVEGetSingleCellDivergence(hcamix,200,10,15);
llshcamix2 = DSAVEGetSingleCellDivergence(hcamix,200,10,15);
figure
scatter(llshcamix1,llshcamix2,3);
xlabel('Log likelihood run 1')
ylabel('Log likelihood run 2')
title('Difference in Cell Divergence Between Two Runs');
pointsToRem = isnan(llshcamix1) | isnan(llshcamix2);
c1 = llshcamix1(~pointsToRem);
c2 = llshcamix2(~pointsToRem);
set(gca,'FontSize',11);

corrcoef(c1,c2)

%% Experiments on why some cells in real datasets have a higher likelihood than the SNO cells
%{
meanAll = mean(hcatSub2.data,2);
hcat500Best = hcatSub2.cellSubset(hcati(1,2001:2500));
mean500Best = mean(hcat500Best.data,2);

[meanAllSort, i500best] = sort(meanAll);
%mean500BestSort = mean500Best(i500best);
mean500BestSort = sort(mean500Best);
xxxx = 1:size(mean500BestSort,1);
figure
semilogy(xxxx,meanAllSort.');
hold on
semilogy(xxxx,mean500BestSort.');
legend({'all','500best'})
%}