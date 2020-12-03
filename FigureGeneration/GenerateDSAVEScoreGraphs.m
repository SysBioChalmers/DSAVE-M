%set random seed to make the results reproducable
rng(1);


%% define datasets

ovasc = SCDep.scd_ovasc;
b10000 = SCDep.scd_pbmcb10000;
bc2 = SCDep.scd_bc2;
[lc,scd_lc_healthy] = SCDep.scd_lc;
pbmc68000 = SCDep.scd_pbmc68000;
hcacb = SCDep.scd_hca_cb;


ovm = ovasc.cellSubset(ovasc.cellType == Celltype.MacrophageOrMonocyte);
bc2t = bc2.cellSubset(bc2.cellType == Celltype.TCellCD4Pos | bc2.cellType == Celltype.TCellCD8Pos | bc2.cellType == Celltype.TCellReg);
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
bc2t_blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_BLOOD'));
lct = lc.cellSubset(lc.cellType == Celltype.TCellCD4Pos | lc.cellType == Celltype.TCellCD8Pos | lc.cellType == Celltype.TCellReg);
lcb = lc.cellSubset(lc.cellType == Celltype.BCell);
lcm = lc.cellSubset(lc.cellType == Celltype.Malignant);
t68000 = pbmc68000.cellSubset(pbmc68000.cellType == Celltype.TCellCD4Pos | pbmc68000.cellType == Celltype.TCellCD8Pos | pbmc68000.cellType == Celltype.TCellReg);
hcat = hcacb.cellSubset(hcacb.cellType == Celltype.TCell);
hcab = hcacb.cellSubset(hcacb.cellType == Celltype.BCell);
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;


lcmix = lct.randSample(3000).innerJoin(lcm.randSample(3000));

%% Fig A

templInfo = DSAVEGetStandardTemplate();
templInfoAllOutliers = templInfo;
templInfoAllOutliers.fractionUpperOutliers = 0;
templInfoAllOutliers.fractionLowerOutliers = 0;


progbar = ProgrBar('DSAVE Score 1: Fig A');


dss = { bc2t.randSample(2500), ...
        ovm, ...
        scd_GSE112845_cd8.randSample(2500) ...
      };
legendsA = { 'BC T cells - mixed pat and tissue','BC T cells - mixed pat and tissue - SNO','OC macrophages','OC macrophages. - SNO','CD8T T cells, blood', 'CD8T T cells, blood - SNO'};
lineStylesA = {'m-','m--','k-','k--','r-','r--'};

numds = size(dss,2);
resdataA = cell(1,numds);

for i = 1:size(dss,2)
    resdataA{1,i} = DSAVECalcBTMScore(dss{1,i}, templInfoAllOutliers, progbar.GetSubContext(0.33), true);
end

figure
for i = 1:numds
    res = resdataA{1,i};
    plot (res.tpms, res.alignedCVs, lineStylesA{1,i*2-1});
    hold on;
    plot (res.tpms, res.samplingCVs, lineStylesA{1,i*2});
    hold on;
end

legend(legendsA)
xlabel('Gene expression (CPM)')
ylabel('ln(CV+1)')
title('Variation per Gene Expression, Unaligned Cell Pop.');
axis([0 1000 0 3]);
set(gca,'FontSize',11);

progbar.Done();

%% Fig B, C and D

progbar = ProgrBar('DSAVE Score 1: Fig B-D');

origdatasets = {ovm, bc2t, bc2t_bc4tumor, bc2t_blood, b10000, scd_GSE112845_cd8};
names = {'OC macr.','OC macr. - SNO','BC T cells - mixed pat and tissue','BC T cells - mixed pat and tissue - SNO','BC tumor T cells, single pat','BC tumor T cells, single pat - SNO','BC blood T cells, single pat','BC blood T cells, single pat - SNO','B10k B cells, blood','B10k B cells, blood - SNO', 'CD8T T cells, blood', 'CD8T T cells, blood - SNO'};

lStyles = {'k','k--','m','m--','b','b--','g','g--','c','c--','r-','r--'};

numds = size(origdatasets,2);
resdata = cell(1,numds);


% Create data

for i = 1:numds
    resdata{1,i} = DSAVECalcBTMScore(origdatasets{1,i}, templInfoAllOutliers, progbar.GetSubContext(1/numds));
end

progbar.Done();

%fig B
legNames = {};
figure
for i = 1:numds
    %skip sample 3 and 4 for this graph
    if i ~= 3 && i ~= 4
        res = resdata{1,i};
        plot (res.tpms, res.alignedCVs, lStyles{1,i*2-1});
        hold on;
        plot (res.tpms, res.samplingCVs, lStyles{1,i*2});
        hold on;
        legNames = [legNames names{1,i*2-1} names{1,i*2}];
    end
end

legend(legNames)
xlabel('Gene expression (CPM)')
ylabel('ln(CV+1)')
title('Variation per Gene Expression, Aligned Cell Pop.');
axis([0 1000 0.8 2.5]);
set(gca,'FontSize',11);

%fig C (zoomed in fig 2)
legNames = {};
figure
for i = 1:numds
    %skip sample 3 and 4 for this graph
    if i ~= 3 && i ~= 4
        res = resdata{1,i};
        plot (res.tpms, res.alignedCVs, lStyles{1,i*2-1});
        hold on;
        plot (res.tpms, res.samplingCVs, lStyles{1,i*2});
        hold on;
        legNames = [legNames names{1,i*2-1} names{1,i*2}];
    end
end

legend(legNames)
xlabel('Gene expression (CPM)')
ylabel('ln(CV+1)')
title('Variation per Gene Expression, Aligned Cell Pop. (Zoomed)');
axis([800 850 0.85 1.3]);
set(gca,'FontSize',11);

%fig D
legNames = {};
figure
for i = 1:numds
    res = resdata{1,i};
    plot (res.tpms, res.differenceCVs, lStyles{1,i*2-1});%use solid styles
    hold on;
    legNames = [legNames names{1,i*2-1}];
end

legend(legNames)
xlabel('Gene expression (CPM)')
ylabel('BTM variation')
title('BTM Variation per Gene Expression');
axis([0 1000 0 0.31]);
set(gca,'FontSize',11);


