
%%Some initial data exploration
figure
histogram(categorical(SCDep.scd_bc2.sampleIds))
[lc,scd_lc_healthy] = SCDep.scd_lc;
figure
histogram(categorical(lc.sampleIds))
title('Whole LC tumor patients');
figure
histogram(categorical(scd_lc_healthy.sampleIds))
title('Whole LC healthy tissue patients');
figure
histogram(categorical(CelltypeId2CelltypeName(lc.paperClass)))
title('Whole LC tumor cell types');

lcpat1 = lc.cellSubset(strcmp(lc.sampleIds,'1'));
figure
histogram(categorical(CelltypeId2CelltypeName(lcpat1.paperClass)))
title('LC pat 1 tumor cell types');


lcpat2 = lc.cellSubset(strcmp(lc.sampleIds,'2'));
figure
histogram(categorical(CelltypeId2CelltypeName(lcpat2.paperClass)))
title('LC pat 2 tumor cell types');

lcpat3 = lc.cellSubset(strcmp(lc.sampleIds,'3'));
figure
histogram(categorical(CelltypeId2CelltypeName(lcpat3.paperClass)))
title('LC pat 3 tumor cell types');

lcpat4 = lc.cellSubset(strcmp(lc.sampleIds,'4'));
figure
histogram(categorical(CelltypeId2CelltypeName(lcpat4.paperClass)))
title('LC pat 4 tumor cell types');

lcpat5 = lc.cellSubset(strcmp(lc.sampleIds,'5'));
figure
histogram(categorical(CelltypeId2CelltypeName(lcpat5.paperClass)))
title('LC pat 5 tumor cell types');

umispat5 = sum(lcpat5.data,1);
mean(umispat5)

lcpat5mal = lcpat5.cellSubset(lcpat5.paperClass == Celltype.Malignant);
mn = TPM(mean(lcpat5mal.data,2));
sum(mn >= 1)
sel5 = mn >= 5;
sum(sel5)
mean(sum(lcpat5mal.data,1))
lcpat5end = lcpat5.cellSubset(lcpat5.paperClass == Celltype.Endothelial);
mnEnd = TPM(mean(lcpat5end.data,2));
mnEndInCanc = mnEnd(sel5,:);
sum(mnEndInCanc > 5)
mean(sum(lcpat5end.data,1))

lcpat5b = lcpat5.cellSubset(lcpat5.paperClass == Celltype.BCell);
mean(sum(lcpat5b.data,1))

healthyEndo = scd_lc_healthy.cellSubset(scd_lc_healthy.paperClass == Celltype.Endothelial);
mean(sum(healthyEndo.data,1))

ovasc = SCDep.scd_ovasc;
bc2 = SCDep.scd_bc2;
pbmc68000 = SCDep.scd_pbmc68000;
ovm = ovasc.cellSubset(ovasc.paperClass == Celltype.MacrophageOrMonocyte);
bc2t = bc2.cellSubset(bc2.paperClass == Celltype.TCellCD4Pos | bc2.paperClass == Celltype.TCellCD8Pos | bc2.paperClass == Celltype.TCellReg);
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
bc2t_blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_BLOOD'));
lct = lc.cellSubset(lc.paperClass == Celltype.TCellCD4Pos | lc.paperClass == Celltype.TCellCD8Pos | lc.paperClass == Celltype.TCellReg);
lcb = lc.cellSubset(lc.paperClass == Celltype.BCell);
lcm = lc.cellSubset(lc.paperClass == Celltype.Malignant);
t68000 = pbmc68000.cellSubset(pbmc68000.paperClass == Celltype.TCellCD4Pos | pbmc68000.paperClass == Celltype.TCellCD8Pos | pbmc68000.paperClass == Celltype.TCellReg);


tdss = {ovm, bc2t, bc2t_bc4tumor, bc2t_blood, b10000, t68000, lct};



%% Fig A - comparison between datasets

%bc2, mixed pat blood t cells
%ovasc macrophages/monocytes, mixed pat
%lc healthy tissue, mixed pat
%livt, mixed pat
%pbmc68000 t cells, pat A
%pbmcb10000
%pbmctcd4mem10000
%hca cb mixed pat t cells
%GSE112845 CD8+ t cells

%from breast cancer
bc2 = SCDep.scd_bc2;
bc2t = bc2.cellSubset(bc2.paperClass == Celltype.TCellCD4Pos | bc2.paperClass == Celltype.TCellCD8Pos | bc2.paperClass == Celltype.TCellReg);
bc2t_bc1blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC1_BLOOD'));
bc2t_bc4blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_BLOOD'));
bc2t_2blood = bc2t_bc1blood.innerJoin(bc2t_bc4blood.randSample(size(bc2t_bc1blood.data,2)));%mix 50% from each patient

%ovarian cancer ascites
ovm = SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.paperClass == Celltype.MacrophageOrMonocyte);

%from lung cancer
[lc,lch] = SCDep.scd_lc;
lcht = lch.cellSubset(lch.paperClass == Celltype.TCellCD4Pos | lch.paperClass == Celltype.TCellCD8Pos | lch.paperClass == Celltype.TCellReg);
lcht_3 = lcht.cellSubset(strcmp(lcht.sampleIds, '3'));
lcht_4 = lcht.cellSubset(strcmp(lcht.sampleIds, '4'));
lcht_5 = lcht.cellSubset(strcmp(lcht.sampleIds, '5'));
num = size(lcht_4.data,2);
lcht_3pat = lcht_4.innerJoin(lcht_3.randSample(num));%mix 33% from each patient
lcht_3pat = lcht_3pat.innerJoin(lcht_5.randSample(num));

%liver cancer
livt = SCDep.scd_livt;%already t cells

%t cells blood
pbmc68000 = SCDep.scd_pbmc68000;
t68000 = pbmc68000.cellSubset(pbmc68000.paperClass == Celltype.TCellCD4Pos | pbmc68000.paperClass == Celltype.TCellCD8Pos | pbmc68000.paperClass == Celltype.TCellReg);

%b cells blood
b10000 = SCDep.scd_pbmcb10000;

%CD4+ memory t cells from blood
pbmctcd4mem10000 = SCDep.scd_pbmctcd4mem10000;

%from human cell atlas cord blood
hcat = hcacb.cellSubset(hcacb.custClass == Celltype.TCell | hcacb.custClass == Celltype.TCellCD4Pos | hcacb.custClass == Celltype.TCellCD8Pos);
hcat_cb1 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB1'));
hcat_cb2 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB2'));
hcat_cb3 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB3'));
num = size(hcat_cb1.data,2);
hcat_3pat = hcat_cb1.innerJoin(hcat_cb2.randSample(num));%33% from each patient
hcat_3pat = hcat_3pat.innerJoin(hcat_cb3.randSample(num));

%GSE112845 CD8+ t cells
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;



dss = { bc2t_2blood, ...
        ovm, ...
        lcht_3pat, ...
        livt, ...
        t68000, ...
        b10000, ...
        pbmctcd4mem10000, ...
        hcat_3pat, ...
        scd_GSE112845_cd8 ...
      };
%names = { 'bc tumor mix 2 pat','bc tumor single pat','bc blood mix 2 pat','bc blood single pat','lc tumor mix 3 pat', 'lc tumor single pat', 'lc ht mix 3 pat', 'lc ht single pat', 'hca cb 3 pat', 'hca cb pat 1', };


numds = size(dss,2);
resdata = cell(1,numds);
scores = zeros(1,numds);

templInfo = DSAVEGetStandardTemplate();

for i = 1:numds
    resdata{1,i} = CalcDSAVE(dss{1,i}, templInfo);
    scores(1,i) = resdata{1,i}.DSAVEScore;
end

disp('Comparison between datasets: Copy into excel sheet');
scores %just copy these values into the excel sheet



%% Fig B - patient to patient variation

%from breast cancer
bc2 = SCDep.scd_bc2;
bc2t = bc2.cellSubset(bc2.paperClass == Celltype.TCellCD4Pos | bc2.paperClass == Celltype.TCellCD8Pos | bc2.paperClass == Celltype.TCellReg);
bc2t_bc1tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC1_TUMOR'));
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
bc2t_2tumors = bc2t_bc1tumor.innerJoin(bc2t_bc4tumor.randSample(size(bc2t_bc1tumor.data,2)));%mix 50% from each patient
bc2t_bc1blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC1_BLOOD'));
bc2t_bc4blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_BLOOD'));
bc2t_2blood = bc2t_bc1blood.innerJoin(bc2t_bc4blood.randSample(size(bc2t_bc1blood.data,2)));%mix 50% from each patient

%from lung cancer
[lc,lch] = SCDep.scd_lc;
lct = lc.cellSubset(lc.paperClass == Celltype.TCellCD4Pos | lc.paperClass == Celltype.TCellCD8Pos | lc.paperClass == Celltype.TCellReg);
lcht = lch.cellSubset(lch.paperClass == Celltype.TCellCD4Pos | lch.paperClass == Celltype.TCellCD8Pos | lch.paperClass == Celltype.TCellReg);
lct_t3 = lct.cellSubset(strcmp(lct.sampleIds, '3'));
lct_t4 = lct.cellSubset(strcmp(lct.sampleIds, '4'));
lct_t5 = lct.cellSubset(strcmp(lct.sampleIds, '5'));
num = size(lct_t3.data,2);
lct_3tumors = lct_t3.innerJoin(lct_t4.randSample(num));%mix 33% from each patient
lct_3tumors = lct_3tumors.innerJoin(lct_t5.randSample(num));
lcht_3 = lcht.cellSubset(strcmp(lcht.sampleIds, '3'));
lcht_4 = lcht.cellSubset(strcmp(lcht.sampleIds, '4'));
lcht_5 = lcht.cellSubset(strcmp(lcht.sampleIds, '5'));
num = size(lcht_4.data,2);
lcht_3pat = lcht_4.innerJoin(lcht_3.randSample(num));%mix 33% from each patient
lcht_3pat = lcht_3pat.innerJoin(lcht_5.randSample(num));

%from human cell atlas cord blood
hcacb = SCDep.scd_hca_cb;
hcat = hcacb.cellSubset(hcacb.custClass == Celltype.TCell | hcacb.custClass == Celltype.TCellCD4Pos | hcacb.custClass == Celltype.TCellCD8Pos);
hcat_cb1 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB1'));
hcat_cb2 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB2'));
hcat_cb3 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB3'));
num = size(hcat_cb1.data,2);
hcat_3pat = hcat_cb1.innerJoin(hcat_cb2.randSample(num));%33% from each patient
hcat_3pat = hcat_3pat.innerJoin(hcat_cb3.randSample(num));

dss = { bc2t_bc1tumor, ...
        bc2t_bc4tumor, ...
        bc2t_2tumors, ...
        lct_t3, ...
        lct_t4, ...
        lct_t5, ...
        lct_3tumors, ...
        hcat_cb1, ...
        hcat_cb2, ...
        hcat_cb3, ...
        hcat_3pat ...
      };

numds = size(dss,2);
resdata = cell(1,numds);
scores = zeros(1,numds);

%templInfo = DSAVEGetStandardTemplate();
%Generate a special template with only 1941 cells to be able to show both
%breast cancer patients
ovm = SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.custClass == Celltype.MacrophageOrMonocyte);
bc2t = SCDep.scd_bc2.cellSubset(SCDep.scd_bc2.paperClass == Celltype.TCellCD4Pos | SCDep.scd_bc2.paperClass == Celltype.TCellCD8Pos | SCDep.scd_bc2.paperClass == Celltype.TCellReg);
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
b10000 = SCDep.scd_pbmcb10000;
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

datasets = {ovm,bc2t_bc4tumor, b10000, scd_GSE112845_cd8};
templInfoSpec = DSAVEGenerateTemplateInfo(bc2t_bc4tumor, datasets, 1941, 570, 0.025, 0.025);

for i = 1:numds
    resdata{1,i} = CalcDSAVE(dss{1,i}, templInfoSpec);
    scores(1,i) = resdata{1,i}.DSAVEScore;
end

disp('Patient to patient: Copy into excel sheet');
scores %just copy these values into the excel sheet

%% Fig C - effect of mixing cell types

%ovasc = SCDep.scd_ovasc;
%b10000 = SCDep.scd_pbmcb10000;
%pbmc68000 = SCDep.scd_pbmc68000;
%hcacb = SCDep.scd_hca_cb;

%add to test SynchronizeGenes
%a = {ovasc, bc2, lc, hcacb, pbmc68000, b10000};


%add to test SynchronizeGenes - sum should be 0
%b = SynchronizeGenes(a,[], true);
%sum(~strcmp(b{1,1}.genes, ovasc.genes))

%from breast cancer
bc2 = SCDep.scd_bc2;
bc2t = bc2.cellSubset(bc2.paperClass == Celltype.TCellCD4Pos | bc2.paperClass == Celltype.TCellCD8Pos | bc2.paperClass == Celltype.TCellReg);
%show cell types in bc2_bc4blood
%figure
%histogram(categorical(CelltypeId2CelltypeName(bc2.cellSubset(strcmp(bc2t.sampleIds, 'BC4_BLOOD')).paperClass)));
bc2_bc4blood = bc2.cellSubset(strcmp(bc2.sampleIds, 'BC4_BLOOD'));
bc2_bc4blood_t = bc2_bc4blood.cellSubset(bc2_bc4blood.paperClass == Celltype.TCellCD4Pos | bc2_bc4blood.paperClass == Celltype.TCellCD8Pos | bc2_bc4blood.paperClass == Celltype.TCellReg);
bc2_bc4blood_b = bc2_bc4blood.cellSubset(bc2_bc4blood.paperClass == Celltype.BCell);
numb = size(bc2_bc4blood_b.data,2);
bc2_bc4blood_mixbc50_50 = bc2_bc4blood_b.innerJoin(bc2_bc4blood_t.randSample(numb));%mix 50% from each patient

%from human cell atlas cord blood
hcacb = SCDep.scd_hca_cb;
hcat = hcacb.cellSubset(hcacb.custClass == Celltype.TCell | hcacb.custClass == Celltype.TCellCD4Pos | hcacb.custClass == Celltype.TCellCD8Pos);
hcatCD8 = hcacb.cellSubset(hcacb.custClass == Celltype.TCellCD8Pos);
hcab = hcacb.cellSubset(hcacb.custClass == Celltype.BCell);
hcat_cb1 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB1'));
hcat_cb2 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB2'));
hcat_cb3 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB3'));
hcab_cb1 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB1'));
hcab_cb2 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB2'));
hcab_cb3 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB3'));
num = size(hcat_cb1.data,2);
numb1 = size(hcab_cb1.data,2);
hcat_1bt50_50 = hcab_cb1.innerJoin(hcat_cb1.randSample(numb1));%50/50 b and t

%from lc
[lc,~] = SCDep.scd_lc;
lcpat5 = lc.cellSubset(strcmp(lc.sampleIds,'5'));
lcp5t = lcpat5.cellSubset(lcpat5.paperClass == Celltype.TCellCD4Pos | lcpat5.paperClass == Celltype.TCellCD8Pos | lcpat5.paperClass == Celltype.TCellReg |  lcpat5.paperClass == Celltype.TCell);
lcp5b = lcpat5.cellSubset(lcpat5.paperClass == Celltype.BCell);
numb = size(lcp5b.data,2);
lc50_50 = lcp5b.innerJoin(lcp5t.randSample(numb));%50/50 b and t


%figure
%histogram(categorical(hcab.sampleIds));


dss = { bc2_bc4blood_t, ...
        bc2_bc4blood_mixbc50_50, ...
        lcp5t, ...
        lcp5b, ...
        lc50_50 ...
        hcat_cb1, ...
        hcab_cb1, ...
        hcat_1bt50_50 ...
      };
%names = { 'bc tumor mix 2 pat','bc tumor single pat','bc blood mix 2 pat','bc blood single pat','lc tumor mix 3 pat', 'lc tumor single pat', 'lc ht mix 3 pat', 'lc ht single pat', 'hca cb 3 pat', 'hca cb pat 1', };


numds = size(dss,2);
resdata = cell(1,numds);
scores = zeros(1,numds);

templInfo = DSAVEGetStandardTemplate();

for i = 1:numds
    resdata{1,i} = CalcDSAVE(dss{1,i}, templInfo);
    scores(1,i) = resdata{1,i}.DSAVEScore;
end

disp('Mixing cell types: Copy into excel sheet');
scores %just copy these values into the excel sheet


%% Supplementary 3A - technical validation using generated cell populations with added noise

ovasc = SCDep.scd_ovasc;
b10000 = SCDep.scd_pbmcb10000;
pbmc68000 = SCDep.scd_pbmc68000;
bc2 = SCDep.scd_bc2;
hcacb = SCDep.scd_hca_cb;

bc2t = bc2.cellSubset(bc2.paperClass == Celltype.TCellCD4Pos | bc2.paperClass == Celltype.TCellCD8Pos | bc2.paperClass == Celltype.TCellReg);
hcat = hcacb.cellSubset(hcacb.custClass == Celltype.TCell | hcacb.custClass == Celltype.TCellCD4Pos | hcacb.custClass == Celltype.TCellCD8Pos);
ovm = ovasc.cellSubset(ovasc.paperClass == Celltype.MacrophageOrMonocyte);
t68000 = pbmc68000.cellSubset(pbmc68000.paperClass == Celltype.TCellCD4Pos | pbmc68000.paperClass == Celltype.TCellCD8Pos | pbmc68000.paperClass == Celltype.TCellReg);


%figure
%histogram(categorical(hcab.sampleIds));

%throw away some cells to save computational time
hcatc = hcat.randSample(10000);
bc2tc = bc2t.randSample(10000);
t68000c = t68000.randSample(10000);

bc2tGenNoise = GenerateSamplingSSDataset(bc2tc, size(bc2tc.data, 2),0.7);
hcatGenNoise = GenerateSamplingSSDataset(hcatc, size(hcatc.data, 2),0.9);
ovmGenNoise = GenerateSamplingSSDataset(ovm, size(ovm.data, 2),1);
t68000GenNoise = GenerateSamplingSSDataset(t68000c, size(t68000c.data, 2),1.5);
b10000GenNoise = GenerateSamplingSSDataset(b10000, size(b10000.data, 2),2);

bc2tGenNoise1 = GenerateSamplingSSDataset(bc2tc, size(bc2tc.data, 2),1);
hcatGenNoise1 = GenerateSamplingSSDataset(hcatc, size(hcatc.data, 2),1);
t68000GenNoise1 = GenerateSamplingSSDataset(t68000c, size(t68000c.data, 2),1);
b10000GenNoise1 = GenerateSamplingSSDataset(b10000, size(b10000.data, 2),1);



dss = { bc2tGenNoise, ...
        hcatGenNoise, ...
        bc2tGenNoise1, ...
        hcatGenNoise1, ...
        ovmGenNoise, ...
        t68000GenNoise1, ...
        b10000GenNoise1, ...
        t68000GenNoise, ...
        b10000GenNoise ...
      };
%names = { 'bc tumor mix 2 pat','bc tumor single pat','bc blood mix 2 pat','bc blood single pat','lc tumor mix 3 pat', 'lc tumor single pat', 'lc ht mix 3 pat', 'lc ht single pat', 'hca cb 3 pat', 'hca cb pat 1', };


numds = size(dss,2);
resdata = cell(1,numds);
scores = zeros(1,numds);

templInfo = DSAVEGetStandardTemplate();

for i = 1:numds
    resdata{1,i} = CalcDSAVE(dss{1,i}, templInfo);
    scores(1,i) = resdata{1,i}.DSAVEScore;
end

disp('Technical validation: Copy into excel sheet');
scores %just copy these values into the excel sheet

%% Fig 3B Supplementary - Mix of monocytes and T cells
hcacb = SCDep.scd_hca_cb;
hca_cb1 = hcacb.cellSubset(strcmp(hcacb.sampleIds,'CB1'));
hcat = hca_cb1.cellSubset(hca_cb1.custClass == Celltype.TCell | hca_cb1.custClass == Celltype.TCellCD4Pos | hca_cb1.custClass == Celltype.TCellCD8Pos);
hcam = hca_cb1.cellSubset(hca_cb1.custClass == Celltype.Monocyte);
hcat2500 = hcat.randSample(2500);
hcam2500 = hcam.randSample(2500);
templInfo = DSAVEGetStandardTemplate();

xes = 0:5:100;
vals = zeros(1,21);

for i = 1:21
    if i == 0
        ds = hcat2500;
    elseif i == 21
        ds = hcam2500;
    else
        frac = (i-1)*0.05;
        ds = hcat2500.randSample(round(2500*frac)).innerJoin(hcam2500.randSample(round(2500*(1-frac))));
    end
    score = CalcDSAVE(ds, templInfo);
    vals(1,i) = score.DSAVEScore;
end

figure
plot(xes,vals);
xlabel('Percent monocytes in mix');
ylabel('DSAVE variation score');
title('Variation per Cell Type Fractions');
set(gca,'FontSize',11);
axis([0 100 0 0.08]);


