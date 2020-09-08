    
%% Fig A and B

templInfo = DSAVEGetStandardTemplate();

progbar = ProgrBar('Cell-wise: Fig A and B');

%breast cancer
bc2 = SCDep.scd_bc2;
bc2t = bc2.cellSubset(bc2.cellType == Celltype.TCellCD4Pos | bc2.cellType == Celltype.TCellCD8Pos | bc2.cellType == Celltype.TCellReg);
bc2t_bc2ln = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC2_LYMPHNODE'));
bc2tSub = bc2t_bc2ln.randSample(2500);
[llsbc2, genesbc2, genellsbc2] = DSAVEGetSingleCellDivergence(bc2tSub, 200, progbar.GetSubContext(0.10));

%plot against number of umis
[bcx,bci] = sort(llsbc2);

%get DSAVE score for the original dataset and for the one where the 500
%most divergent have been removed
bcDSAVEScore = DSAVECalcBTMScore(bc2tSub, templInfo, progbar.GetSubContext(0.06));

%now remove the 500 worst cells
bc2tSub500Less = bc2tSub.cellSubset(bci(1,501:2500));
bcDSAVEScore500Less = DSAVECalcBTMScore(bc2tSub500Less, templInfo, progbar.GetSubContext(0.06));

%now for the SNO
bcSNO = DSAVEGenerateSNODataset(bc2tSub, progbar.GetSubContext(0.01));
[llsbcSNO, genesbc2SNO, genellsbc2SNO] = DSAVEGetSingleCellDivergence(bcSNO, 200, progbar.GetSubContext(0.10));
[bcsnox,bcsnoi] = sort(llsbcSNO);

%hca
hcacb = SCDep.scd_hca_cb;
hca_cb1 = hcacb.cellSubset(strcmp(hcacb.sampleIds,'CB1'));
hcat = hca_cb1.cellSubset(hca_cb1.cellType == Celltype.TCell | hca_cb1.cellType == Celltype.TCellCD4Pos | hca_cb1.cellType == Celltype.TCellCD8Pos);
hcab = hca_cb1.cellSubset(hca_cb1.cellType == Celltype.BCell);
hcam = hca_cb1.cellSubset(hca_cb1.cellType == Celltype.Monocyte);
hcatSub2 = hcat.randSample(2500);


llshcat = DSAVEGetSingleCellDivergence(hcatSub2, 200, progbar.GetSubContext(0.10));
hcatDSAVEScore = DSAVECalcBTMScore(hcatSub2, templInfo, progbar.GetSubContext(0.06));
[hcatx,hcati] = sort(llshcat);
hcat500Less = hcatSub2.cellSubset(hcati(1,501:2500));
hcat500LessDSAVEScore = DSAVECalcBTMScore(hcat500Less, templInfo, progbar.GetSubContext(0.06));

hcaSNO = DSAVEGenerateSNODataset(hcatSub2, progbar.GetSubContext(0.01));
llshcaSNO = DSAVEGetSingleCellDivergence(hcaSNO, 200, progbar.GetSubContext(0.10));
[hcatsnox,hcatsnoi] = sort(llshcaSNO);


%B10k
b10k = SCDep.scd_pbmcb10000;
b10kSub = b10k.randSample(2500);
llsb10k = DSAVEGetSingleCellDivergence(b10kSub, 200, progbar.GetSubContext(0.10));
[b10kx,b10ki] = sort(llsb10k);
%sno
b10kSNO = DSAVEGenerateSNODataset(b10kSub, progbar.GetSubContext(0.01));
llsb10kSNO = DSAVEGetSingleCellDivergence(b10kSNO, 200, progbar.GetSubContext(0.10));
[b10ksnox,b10ksnoi] = sort(llsb10kSNO);

%get DSAVE score for the original dataset and for the one where the 500
%most divergent have been removed
b10kDSAVEScore = DSAVECalcBTMScore(b10kSub, templInfo, progbar.GetSubContext(0.06));

%now remove the 500 worst cells
b10k500Less = b10kSub.cellSubset(b10ki(1,501:2500));
b10kDSAVEScore500Less = DSAVECalcBTMScore(b10k500Less, templInfo, progbar.GetSubContext(0.06));


%Now plot fig A:
xs = 1:2500;
figure
plot(xs, -bcx, '-', 'Color', [0, 0.4470, 0.7410],'LineWidth',2);
hold on
plot(xs, -bcsnox, '--', 'Color', [0, 0.4470, 0.7410],'LineWidth',2);
hold on
plot(xs, -b10kx, '-', 'Color', [0.8500, 0.3250, 0.0980],'LineWidth',2);
hold on
plot(xs, -b10ksnox, '--', 'Color', [0.8500, 0.3250, 0.0980],'LineWidth',2);
hold on
plot(xs, -hcatx, '-', 'Color',[0.800, 0.600, 0.100],'LineWidth',2);
hold on
plot(xs, -hcatsnox, '--', 'Color',[0.800, 0.600, 0.100],'LineWidth',2);
hold on
xlabel('Cell index')
ylabel('Cell divergence')
title('Cell Divergence');
h = legend({'BC LN T cells, single pat', 'BC LN T cells, single pat - SNO', 'B10k B cells, single pat', 'B10k B cells, single pat - SNO', 'HCA CB T cells, single pat', 'HCA CB T cells, single pat - SNO'});
set(h, 'Position', [0.4637    0.647    0.4179    0.2583])
set(gca,'FontSize',11);
axis([0, 2500, 500, 3100]);

progbar.Done();

%check if we can explain why the SNO is lower in the right part
% bcx
% genesPerCell = sum(bc2tSub.data > 0, 1);
% genesPerCellSort = genesPerCell(bci);
% figure
% scatter(genesPerCellSort,bcx);
% meanExprTest = TPM(sum(bc2tSub.data,2));
% lowGenes = meanExprTest < 10;
% highGenes = meanExprTest > 1000;
% lowGenesPerCell = sum(bc2tSub.data(lowGenes,:) > 0, 1);
% fracLowGenes = lowGenesPerCell ./ genesPerCell;
% fracLowGenesSort = fracLowGenes(bci);
% highGenesPerCell = sum(bc2tSub.data(highGenes,:) > 0, 1);
% fracHighGenes = highGenesPerCell ./ genesPerCell;
% fracHighGenesSort = fracHighGenes(bci);
% figure
% scatter(fracLowGenesSort,bcx);
% lsline
% %compare with SNO
% genesPerCellSno = sum(bcSNO.data > 0, 1);
% lowGenesPerCellSno = sum(bcSNO.data(lowGenes,:) > 0, 1);
% fracLowGenesSno = lowGenesPerCellSno ./ genesPerCellSno;
% fracLowGenesSortSno = fracLowGenesSno(bcsnoi);
% figure
% scatter(fracLowGenesSortSno,bcsnox);
% lsline
%conclusion - it cannot be explained by that some cells in the real datasets has less lowly
%expressed genes - the reason is likely more complex

%Data for Fig B:
disp('Fig 5B data, copy to excel:');
a = [bcDSAVEScore.DSAVEScore bcDSAVEScore500Less.DSAVEScore b10kDSAVEScore.DSAVEScore b10kDSAVEScore500Less.DSAVEScore hcatDSAVEScore.DSAVEScore hcat500LessDSAVEScore.DSAVEScore]


%% Fig C - Num UMIs vs log likelihood
numUMIshcatSub2 = sum(hcatSub2.data,1);

linevalYs = 1500:500:10500;
linevalXes = zeros(1,size(linevalYs,2));
%generate line
for i = 1:size(linevalYs,2)
   ii = linevalYs(i);
   lb = ii - 250;
   ub = ii + 250;
   sel = numUMIshcatSub2 >= lb & numUMIshcatSub2 <= ub;
   linevalXes(1,i) = mean(llshcat(sel));
end

figure
scatter(-llshcat,numUMIshcatSub2,3);
hold on
plot(-linevalXes,linevalYs,'Color',[0, 0.8, 0.0],'LineWidth',2);
xlabel('Cell divergence')
ylabel('UMI counts')
title('UMI Counts vs Cell Divergence');
legend({'Individual cell', 'Mean divergence'});
set(gca,'FontSize',11);
axis([1000 2700 0 14000]);

%% Fig D. MT-genes

progbar = ProgrBar('Cell-wise: Fig D');

matches = regexp(bc2tSub.genes, '^(MT-.*)$', 'match');
a = cellfun(@isempty, matches);
b = ~a;
sum(b)
bc2tSub.genes(b)

sumMito = sum(bc2tSub.data(b,:),1);
sumUMIs = sum(bc2tSub.data,1);
percMito = sumMito ./ sumUMIs;

bc2tSubNoMT = bc2tSub.geneSubset(a);

%calc divergence for the dataset when the MT genes have been removed to
%make sure that the divergence is not affected by the MT genes in
%themselves.
llsbc2NoMT = DSAVEGetSingleCellDivergence(bc2tSubNoMT, 200, progbar.GetSubContext(1));
progbar.Done();


linevalYs = 0.01:0.01:0.2;
linevalXes = zeros(1,20);
%generate line
for i = 1:20;
   ii = i*.01;
   lb = ii - 0.01;
   ub = ii + 0.01;
   sel = percMito >= lb & percMito <= ub;
   linevalXes(1,i) = mean(llsbc2NoMT(sel));
end


figure
scatter(-llsbc2,percMito,3);
hold on
plot(-linevalXes,linevalYs,'Color',[0, 0.8, 0.0],'LineWidth',2);
xlabel('Cell divergence')
ylabel('Fraction mitochondrial gene counts')
title('Fraction MT Counts vs Cell Divergence');
legend({'Individual cell', 'Mean divergence'});
set(gca,'FontSize',11);
axis([800 1600 0 0.25]);

%% Fig 5B in supplementary - Number of detected genes
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
scatter(-llshcat,detGenes,3);
hold on
plot(-linevalXes,linevalYs,'Color',[0, 0.8, 0.0],'LineWidth',2);
xlabel('Cell divergence')
ylabel('Number of detected genes')
title('Number of Detected Genes vs Cell Divergence');
legend({'Individual cell', 'Mean divergence'});
set(gca,'FontSize',11);
axis([1000 2700 450 2800]);


%% Fig 5A in supplementary. Technical verification - run1 vs run 2.
%run divergence 2 times and plot
progbar = ProgrBar('Cell-wise: Fig 4A Suppl.');

hcat2500 = hcat.randSample(2500);
hcam2500 = hcam.randSample(2500);

hcatSub = hcat2500.randSample(2250);
hcamSub = hcam2500.randSample(250);
hcamix = hcatSub.innerJoin(hcamSub);%the 500 first cells will be b cells


llshcamix1 = DSAVEGetSingleCellDivergence(hcamix,200,progbar.GetSubContext(0.48),10,15);
llshcamix2 = DSAVEGetSingleCellDivergence(hcamix,200,progbar.GetSubContext(0.48),10,15);
figure
scatter(-llshcamix1,-llshcamix2,3);
xlabel('Divergence run 1')
ylabel('Divergence run 2')
title('Difference in Cell Divergence Between Two Runs');
pointsToRem = isnan(llshcamix1) | isnan(llshcamix2);
c1 = llshcamix1(~pointsToRem);
c2 = llshcamix2(~pointsToRem);
set(gca,'FontSize',11);

progbar.Done();

disp('Correlation coefficient:');
corrcoef(c1,c2)


