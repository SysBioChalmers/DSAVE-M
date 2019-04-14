

%% define datasets

ovasc = SCDep.scd_ovasc;
b10000 = SCDep.scd_pbmcb10000;
bc2 = SCDep.scd_bc2;
[lc,scd_lc_healthy] = SCDep.scd_lc;
pbmc68000 = SCDep.scd_pbmc68000;
hcacb = SCDep.scd_hca_cb;

%add to test SynchronizeGenes
%a = {ovasc, bc2, lc, hcacb, pbmc68000, b10000};


%add to test SynchronizeGenes - sum should be 0
%b = SynchronizeGenes(a,[], true);
%sum(~strcmp(b{1,1}.genes, ovasc.genes))


ovm = ovasc.cellSubset(ovasc.paperClass == Celltype.MacrophageOrMonocyte);
bc2t = bc2.cellSubset(bc2.paperClass == Celltype.TCellCD4Pos | bc2.paperClass == Celltype.TCellCD8Pos | bc2.paperClass == Celltype.TCellReg);
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
bc2t_blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_BLOOD'));
lct = lc.cellSubset(lc.paperClass == Celltype.TCellCD4Pos | lc.paperClass == Celltype.TCellCD8Pos | lc.paperClass == Celltype.TCellReg);
lcb = lc.cellSubset(lc.paperClass == Celltype.BCell);
lcm = lc.cellSubset(lc.paperClass == Celltype.Malignant);
t68000 = pbmc68000.cellSubset(pbmc68000.paperClass == Celltype.TCellCD4Pos | pbmc68000.paperClass == Celltype.TCellCD8Pos | pbmc68000.paperClass == Celltype.TCellReg);
hcat = hcacb.cellSubset(hcacb.custClass == Celltype.TCell);
hcab = hcacb.cellSubset(hcacb.custClass == Celltype.BCell);
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;


lcmix = lct.randSample(3000).innerJoin(lcm.randSample(3000));

%% Fig A

ub = 100000;
lb = 0.5;
n = 6000;

dss = { bc2t.randSample(2500), ...
        ovm, ...
        scd_GSE112845_cd8.randSample(2500) ...
      };
legendsA = { 'BC tumor T cells','BC tumor T cells - SNO','OC macr.','OC macr. - SNO','CD8+ T cells, blood', 'CD8+ T cells, blood - SNO'};
lineStylesA = {'m-','m--','k-','k--','r-','r--'};

numds = size(dss,2);
resdataA = cell(1,numds);

templInfo = DSAVEGetStandardTemplate();
templInfoAllOutliers = templInfo;
templInfoAllOutliers.fractionUpperOutliers = 0;
templInfoAllOutliers.fractionLowerOutliers = 0;

for i = 1:size(dss,2)
    resdataA{1,i} = CalcDSAVE(dss{1,i}, templInfoAllOutliers, true);
end

%legNames = {};
figure
for i = 1:numds
    res = resdataA{1,i};
    plot (res.tpms, res.alignedCVs, lineStyles{1,i*2-1});
    hold on;
    plot (res.tpms, res.samplingCVs, lineStyles{1,i*2});
    hold on;
%    legNames = [legNames names{1,i*3-2} names{1,i*3-1}];
end

legend(legends)
xlabel('Gene expression (CPM)')
ylabel('Log_2(CV+1)')
title('Variation per Gene Expression, Unaligned Cell Pop.');
axis([0 1000 0 4]);
set(gca,'FontSize',11);

%{
figure
%for some reason it doesn't work with the linestyles using a vector,
%so we'll have to loop
for i = 1:size(X2,2)
    h = plot(X2(:,i), Y2(:,i), lineStyles{1,i});
    hold on
end
%h = plot(X, Y, 'linestyle',{'--','--','--','--'});
%hbc = get(h, 'Children');
%set(hbc{1}, 'FaceColor', 'r');
legend(legends);
xlabel('Pool size (number of cells)')
ylabel('Log_2(CV+1)')
title('Variation per Cell Pool Size, Non-aligned Cell Pop.');
axis([0 6000 0 1.2]);
set(gca,'FontSize',11);
%}

%% Fig B, C and D

origdatasets = {ovm, bc2t, bc2t_bc4tumor, bc2t_blood, b10000, scd_GSE112845_cd8};
names = {'OC macr.','OC macr. - SNO','OC macr.','BC tumor T cells, mixed pat','BC tumor T cells, mixed pat - SNO','BC tumor T cells, mixed pat','BC tumor T cells, single pat','BC tumor T cells, single pat - SNO','BC tumor T cells, single pat','BC blood T cells, mixed pat','BC blood T cells, mixed pat - SNO','BC blood T cells, mixed pat','B cells, blood','B cells, blood - SNO','B cells, blood', 'CD8+ T cells, blood', 'CD8+ T cells, blood - SNO', 'CD8+ T cells, blood'};

lStyles = {'k','k--','k:','m','m--','m:','b','b--','b:','g','g--','g:','c','c--','c:','r-','r--','r:'};

numds = size(origdatasets,2);
resdata = cell(1,numds);
scores = zeros(1,numds);


% Create data

for i = 1:numds
    resdata{1,i} = CalcDSAVE(origdatasets{1,i}, templInfoAllOutliers);
    scores(1,i) = resdata{1,i}.DSAVEScore;
end



%fig B
legNames = {};
figure
for i = 1:numds
    %skip sample 3 and 4 for this graph
    if i ~= 3 & i ~= 4
        res = resdata{1,i};
        plot (res.tpms, res.alignedCVs, lStyles{1,i*3-2});
        hold on;
        plot (res.tpms, res.samplingCVs, lStyles{1,i*3-1});
        hold on;
        legNames = [legNames names{1,i*3-2} names{1,i*3-1}];
    end
end

legend(legNames)
xlabel('Gene expression (CPM)')
ylabel('Log_2(CV+1)')
title('Variation per Gene Expression, Aligned Cell Pop.');
axis([0 1000 1 3]);
set(gca,'FontSize',11);

%fig C (zoomed in fig 2)
legNames = {};
figure
for i = 1:numds
    %skip sample 3 and 4 for this graph
    if i ~= 3 & i ~= 4
        res = resdata{1,i};
        plot (res.tpms, res.alignedCVs, lStyles{1,i*3-2});
        hold on;
        plot (res.tpms, res.samplingCVs, lStyles{1,i*3-1});
        hold on;
        legNames = [legNames names{1,i*3-2} names{1,i*3-1}];
    end
end

legend(legNames)
xlabel('Gene expression (CPM)')
ylabel('Log_2(CV+1)')
title('Variation per Gene Expression, Aligned Cell Pop. (Zoomed)');
axis([400 450 1.5 1.9]);
set(gca,'FontSize',11);

%fig D
legNames = {};
figure
for i = 1:numds
    res = resdata{1,i};
    plot (res.tpms, res.differenceCVs, lStyles{1,i*3-2});%use solid styles
    hold on;
    legNames = [legNames names{1,i*3}];
end

legend(legNames)
xlabel('Gene expression (CPM)')
ylabel('BTM variation')
title('BTM Variation per Gene Expression');
axis([0 1000 0 0.41]);
set(gca,'FontSize',11);

