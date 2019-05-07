
%% LC
[lc,scd_lc_healthy] = SCDep.scd_lc;
lcpat1 = lc.cellSubset(strcmp(lc.sampleIds,'1'));
lcpat2 = lc.cellSubset(strcmp(lc.sampleIds,'2'));
lcpat3 = lc.cellSubset(strcmp(lc.sampleIds,'3'));
lcpat4 = lc.cellSubset(strcmp(lc.sampleIds,'4'));
lcpat5 = lc.cellSubset(strcmp(lc.sampleIds,'5'));
lchpat3 = scd_lc_healthy.cellSubset(strcmp(scd_lc_healthy.sampleIds,'3'));
lchpat4 = scd_lc_healthy.cellSubset(strcmp(scd_lc_healthy.sampleIds,'4'));
lchpat5 = scd_lc_healthy.cellSubset(strcmp(scd_lc_healthy.sampleIds,'5'));


%{
figure
histogram(categorical(lc.sampleIds))
title('Whole LC tumor patients');
figure
histogram(categorical(scd_lc_healthy.sampleIds))
title('Whole LC healthy tissue patients');

figure
histogram(categorical(CelltypeId2CelltypeName(lc.cellType)))
title('Whole LC tumor cell types');

figure
histogram(categorical(CelltypeId2CelltypeName(lcpat1.cellType)))
title('LC pat 1 tumor cell types');

figure
histogram(categorical(CelltypeId2CelltypeName(lcpat2.cellType)))
title('LC pat 2 tumor cell types');

figure
histogram(categorical(CelltypeId2CelltypeName(lcpat3.cellType)))
title('LC pat 3 tumor cell types');

figure
histogram(categorical(CelltypeId2CelltypeName(lcpat4.cellType)))
title('LC pat 4 tumor cell types');

figure
histogram(categorical(CelltypeId2CelltypeName(lcpat5.cellType)))
title('LC pat 5 tumor cell types');

figure
histogram(categorical(CelltypeId2CelltypeName(lchpat3.cellType)))
title('LC healthy pat 3 tumor cell types');

figure
histogram(categorical(CelltypeId2CelltypeName(lchpat4.cellType)))
title('LC healthy pat 4 tumor cell types');

figure
histogram(categorical(CelltypeId2CelltypeName(lchpat5.cellType)))
title('LC healthy pat 5 tumor cell types');

umispat5 = sum(lcpat5.data,1);
mean(umispat5)

lcpat5mal = lcpat5.cellSubset(lcpat5.cellType == Celltype.Malignant);
mn = TPM(mean(lcpat5mal.data,2));
sum(mn >= 1)
sel5 = mn >= 5;
sum(sel5)
mean(sum(lcpat5mal.data,1))
lcpat5end = lcpat5.cellSubset(lcpat5.cellType == Celltype.Endothelial);
mnEnd = TPM(mean(lcpat5end.data,2));
mnEndInCanc = mnEnd(sel5,:);
sum(mnEndInCanc > 5)
mean(sum(lcpat5end.data,1))

lcpat5b = lcpat5.cellSubset(lcpat5.cellType == Celltype.BCell);
mean(sum(lcpat5b.data,1))

healthyEndo = scd_lc_healthy.cellSubset(scd_lc_healthy.cellType == Celltype.Endothelial);
mean(sum(healthyEndo.data,1))
%}

%create all LC subpopulations
lchp5t = lchpat5.cellSubset(lchpat5.cellType == Celltype.TCellCD4Pos | lchpat5.cellType == Celltype.TCellCD8Pos | lchpat5.cellType == Celltype.TCellReg |  lchpat5.cellType == Celltype.TCell );
lchp5t.name = 'lchp5t';
lchp5m = lchpat5.cellSubset(lchpat5.cellType == Celltype.Macrophage);
lchp5m.name = 'lchp5mac';
lchp4t = lchpat4.cellSubset(lchpat4.cellType == Celltype.TCellCD4Pos | lchpat4.cellType == Celltype.TCellCD8Pos | lchpat4.cellType == Celltype.TCellReg |  lchpat4.cellType == Celltype.TCell );
lchp4t.name = 'lchp4t';
lchp3t = lchpat3.cellSubset(lchpat3.cellType == Celltype.TCellCD4Pos | lchpat3.cellType == Celltype.TCellCD8Pos | lchpat3.cellType == Celltype.TCellReg |  lchpat3.cellType == Celltype.TCell );
lchp3t.name = 'lchp3t';
lcp5t = lcpat5.cellSubset(lcpat5.cellType == Celltype.TCellCD4Pos | lcpat5.cellType == Celltype.TCellCD8Pos | lcpat5.cellType == Celltype.TCellReg |  lcpat5.cellType == Celltype.TCell);
lcp5t.name = 'lcp5t';
lcp5mal = lcpat5.cellSubset(lcpat5.cellType == Celltype.Malignant);
lcp5mal.name = 'lcp5mal';
lcp5b = lcpat5.cellSubset(lcpat5.cellType == Celltype.BCell);
lcp5b.name = 'lcp5b';
lcp3t = lcpat3.cellSubset(lcpat3.cellType == Celltype.TCellCD4Pos | lcpat3.cellType == Celltype.TCellCD8Pos | lcpat3.cellType == Celltype.TCellReg |  lcpat3.cellType == Celltype.TCell);
lcp3t.name = 'lcp3t';
lcp3mal = lcpat3.cellSubset(lcpat3.cellType == Celltype.Malignant);
lcp3mal.name = 'lcp3mal';
lcp3m = lcpat3.cellSubset(lcpat3.cellType == Celltype.Macrophage);
lcp3m.name = 'lcp3mac';

%%BC
bc = SCDep.scd_bc2;
bc1b = bc.cellSubset(strcmp(bc.sampleIds,'BC1_BLOOD'));
bc1n = bc.cellSubset(strcmp(bc.sampleIds,'BC1_NORMAL'));
bc1t = bc.cellSubset(strcmp(bc.sampleIds,'BC1_TUMOR'));
bc2ln = bc.cellSubset(strcmp(bc.sampleIds,'BC2_LYMPHNODE'));
bc2n = bc.cellSubset(strcmp(bc.sampleIds,'BC2_NORMAL'));
bc2t = bc.cellSubset(strcmp(bc.sampleIds,'BC2_TUMOR'));
bc4b = bc.cellSubset(strcmp(bc.sampleIds,'BC4_BLOOD'));
bc4t = bc.cellSubset(strcmp(bc.sampleIds,'BC4_TUMOR'));
bc5t = bc.cellSubset(strcmp(bc.sampleIds,'BC5_TUMOR'));
bc6t = bc.cellSubset(strcmp(bc.sampleIds,'BC6_TUMOR'));
bc7t = bc.cellSubset(strcmp(bc.sampleIds,'BC7_TUMOR'));
bc8t = bc.cellSubset(strcmp(bc.sampleIds,'BC8_TUMOR'));
%{
figure
histogram(categorical(bc.sampleIds))
figure
histogram(categorical(CelltypeId2CelltypeName(bc1b.cellType)))
title('BC pat 1 blood cell types');
figure
histogram(categorical(CelltypeId2CelltypeName(bc1n.cellType)))
title('BC pat 1 normal cell types');
figure
histogram(categorical(CelltypeId2CelltypeName(bc1t.cellType)))
title('BC pat 1 tumor cell types');
figure
histogram(categorical(CelltypeId2CelltypeName(bc2ln.cellType)))
title('BC pat 2 ln cell types');
figure
histogram(categorical(CelltypeId2CelltypeName(bc2n.cellType)))
title('BC pat 2 normal cell types');%there are no populations with enough cells here
figure
histogram(categorical(CelltypeId2CelltypeName(bc2t.cellType)))
title('BC pat 2 tumor cell types');
figure
histogram(categorical(CelltypeId2CelltypeName(bc4b.cellType)))
title('BC pat 4 blood cell types');
figure
histogram(categorical(CelltypeId2CelltypeName(bc4t.cellType)))
title('BC pat 4 tumor cell types');
figure
histogram(categorical(CelltypeId2CelltypeName(bc5t.cellType)))
title('BC pat 5 tumor cell types');%there are no populations with enough cells here
figure
histogram(categorical(CelltypeId2CelltypeName(bc6t.cellType)))
title('BC pat 6 tumor cell types');
figure
histogram(categorical(CelltypeId2CelltypeName(bc7t.cellType)))
title('BC pat 7 tumor cell types');%there are no populations with enough cells here
figure
histogram(categorical(CelltypeId2CelltypeName(bc8t.cellType)))
title('BC pat 8 tumor cell types');%there are no populations with enough cells here
%}
bc1b_t = bc1b.cellSubset(bc1b.cellType == Celltype.TCellCD4Pos | bc1b.cellType == Celltype.TCellCD8Pos | bc1b.cellType == Celltype.TCellReg |  bc1b.cellType == Celltype.TCell);
bc1b_t.name = 'bc1b_t';
bc1n_t = bc1n.cellSubset(bc1n.cellType == Celltype.TCellCD4Pos | bc1n.cellType == Celltype.TCellCD8Pos | bc1n.cellType == Celltype.TCellReg |  bc1n.cellType == Celltype.TCell);
bc1n_t.name = 'bc1n_t';
bc1t_t = bc1t.cellSubset(bc1t.cellType == Celltype.TCellCD4Pos | bc1t.cellType == Celltype.TCellCD8Pos | bc1t.cellType == Celltype.TCellReg |  bc1t.cellType == Celltype.TCell);
bc1t_t.name = 'bc1t_t';
bc2ln_t = bc2ln.cellSubset(bc2ln.cellType == Celltype.TCellCD4Pos | bc2ln.cellType == Celltype.TCellCD8Pos | bc2ln.cellType == Celltype.TCellReg |  bc2ln.cellType == Celltype.TCell);
bc2ln_t.name = 'bc2ln_t';
bc2ln_b = bc2ln.cellSubset(bc2ln.cellType == Celltype.BCell);
bc2ln_b.name = 'bc2ln_b';
bc2t_t = bc2t.cellSubset(bc2t.cellType == Celltype.TCellCD4Pos | bc2t.cellType == Celltype.TCellCD8Pos | bc2t.cellType == Celltype.TCellReg |  bc2t.cellType == Celltype.TCell);
bc2t_t.name = 'bc2t_t';
bc4b_t = bc4b.cellSubset(bc4b.cellType == Celltype.TCellCD4Pos | bc4b.cellType == Celltype.TCellCD8Pos | bc4b.cellType == Celltype.TCellReg |  bc4b.cellType == Celltype.TCell);
bc4b_t.name = 'bc4b_t';
bc4b_b = bc4b.cellSubset(bc4b.cellType == Celltype.BCell);
bc4b_b.name = 'bc4b_b';
bc4b_nk = bc4b.cellSubset(bc4b.cellType == Celltype.NKCell);
bc4b_nk.name = 'bc4b_nk';
bc4t_t = bc4t.cellSubset(bc4t.cellType == Celltype.TCellCD4Pos | bc4t.cellType == Celltype.TCellCD8Pos | bc4t.cellType == Celltype.TCellReg |  bc4t.cellType == Celltype.TCell);
bc4t_t.name = 'bc4t_t';
bc6t_t = bc6t.cellSubset(bc6t.cellType == Celltype.TCellCD4Pos | bc6t.cellType == Celltype.TCellCD8Pos | bc6t.cellType == Celltype.TCellReg |  bc6t.cellType == Celltype.TCell);
bc6t_t.name = 'bc6t_t';
bc6t_m = bc6t.cellSubset(bc6t.cellType == Celltype.Macrophage);
bc6t_m.name = 'bc6t_m';
bc8t_t = bc8t.cellSubset(bc8t.cellType == Celltype.TCellCD4Pos | bc8t.cellType == Celltype.TCellCD8Pos | bc8t.cellType == Celltype.TCellReg |  bc8t.cellType == Celltype.TCell);
bc8t_t.name = 'bc8t_t';

%% HCA

hca = SCDep.scd_hca_cb;
hca1 = hca.cellSubset(strcmp(hca.sampleIds,'CB1'));
hca2 = hca.cellSubset(strcmp(hca.sampleIds,'CB2'));
hca3 = hca.cellSubset(strcmp(hca.sampleIds,'CB3'));
hca4 = hca.cellSubset(strcmp(hca.sampleIds,'CB4'));
hca5 = hca.cellSubset(strcmp(hca.sampleIds,'CB5'));
hca6 = hca.cellSubset(strcmp(hca.sampleIds,'CB6'));
hca7 = hca.cellSubset(strcmp(hca.sampleIds,'CB7'));
hca8 = hca.cellSubset(strcmp(hca.sampleIds,'CB8'));
%{
figure
histogram(categorical(hca.sampleIds))
figure
histogram(categorical(CelltypeId2CelltypeName(hca1.cellType)))
title('HCA pat 1 cell types');%there are no populations with enough cells here
figure
histogram(categorical(CelltypeId2CelltypeName(hca2.cellType)))
title('HCA pat 2 cell types');%there are no populations with enough cells here
figure
histogram(categorical(CelltypeId2CelltypeName(hca3.cellType)))
title('HCA pat 3 cell types');%there are no populations with enough cells here
figure
histogram(categorical(CelltypeId2CelltypeName(hca4.cellType)))
title('HCA pat 4 cell types');%there are no populations with enough cells here
figure
histogram(categorical(CelltypeId2CelltypeName(hca5.cellType)))
title('HCA pat 5 cell types');%there are no populations with enough cells here
figure
histogram(categorical(CelltypeId2CelltypeName(hca6.cellType)))
title('HCA pat 6 cell types');%there are no populations with enough cells here
figure
histogram(categorical(CelltypeId2CelltypeName(hca7.cellType)))
title('HCA pat 7 cell types');%there are no populations with enough cells here
figure
histogram(categorical(CelltypeId2CelltypeName(hca8.cellType)))
title('HCA pat 8 cell types');%there are no populations with enough cells here
%}
hca1_t = hca1.cellSubset(hca1.cellType == Celltype.TCellCD4Pos | hca1.cellType == Celltype.TCellCD8Pos | hca1.cellType == Celltype.TCellReg |  hca1.cellType == Celltype.TCell);
hca1_t.name = 'hca1_t';
hca1_b = hca1.cellSubset(hca1.cellType == Celltype.BCell);
hca1_b.name = 'hca1_b';
hca1_mon = hca1.cellSubset(hca1.cellType == Celltype.Monocyte);
hca1_mon.name = 'hca1_mon';
hca1_nk = hca1.cellSubset(hca1.cellType == Celltype.NKCell);
hca1_nk.name = 'hca1_nk';
hca2_t = hca2.cellSubset(hca2.cellType == Celltype.TCellCD4Pos | hca2.cellType == Celltype.TCellCD8Pos | hca2.cellType == Celltype.TCellReg |  hca2.cellType == Celltype.TCell);
hca2_t.name = 'hca2_t';
hca2_b = hca2.cellSubset(hca2.cellType == Celltype.BCell);
hca2_b.name = 'hca2_b';
hca2_mon = hca2.cellSubset(hca2.cellType == Celltype.Monocyte);
hca2_mon.name = 'hca2_mon';
hca2_nk = hca2.cellSubset(hca2.cellType == Celltype.NKCell);
hca2_nk.name = 'hca2_nk';
hca3_t = hca3.cellSubset(hca3.cellType == Celltype.TCellCD4Pos | hca3.cellType == Celltype.TCellCD8Pos | hca3.cellType == Celltype.TCellReg |  hca3.cellType == Celltype.TCell);
hca3_t.name = 'hca3_t';
hca3_b = hca3.cellSubset(hca3.cellType == Celltype.BCell);
hca3_b.name = 'hca3_b';
hca4_t = hca4.cellSubset(hca4.cellType == Celltype.TCellCD4Pos | hca4.cellType == Celltype.TCellCD8Pos | hca4.cellType == Celltype.TCellReg |  hca4.cellType == Celltype.TCell);
hca4_t.name = 'hca4_t';
hca4_b = hca4.cellSubset(hca4.cellType == Celltype.BCell);
hca4_b.name = 'hca4_b';
hca4_mon = hca4.cellSubset(hca4.cellType == Celltype.Monocyte);
hca4_mon.name = 'hca4_mon';
hca4_nk = hca4.cellSubset(hca4.cellType == Celltype.NKCell);
hca4_nk.name = 'hca4_nk';
hca5_t = hca5.cellSubset(hca5.cellType == Celltype.TCellCD4Pos | hca5.cellType == Celltype.TCellCD8Pos | hca5.cellType == Celltype.TCellReg |  hca5.cellType == Celltype.TCell);
hca5_t.name = 'hca5_t';
hca5_b = hca5.cellSubset(hca5.cellType == Celltype.BCell);
hca5_b.name = 'hca5_b';
hca5_mon = hca5.cellSubset(hca5.cellType == Celltype.Monocyte);
hca5_mon.name = 'hca5_mon';
hca5_nk = hca5.cellSubset(hca5.cellType == Celltype.NKCell);
hca5_nk.name = 'hca5_nk';
hca6_t = hca6.cellSubset(hca6.cellType == Celltype.TCellCD4Pos | hca6.cellType == Celltype.TCellCD8Pos | hca6.cellType == Celltype.TCellReg |  hca6.cellType == Celltype.TCell);
hca6_t.name = 'hca6_t';
hca6_b = hca6.cellSubset(hca6.cellType == Celltype.BCell);
hca6_b.name = 'hca6_b';
hca6_mon = hca6.cellSubset(hca6.cellType == Celltype.Monocyte);
hca6_mon.name = 'hca6_mon';
hca7_t = hca7.cellSubset(hca7.cellType == Celltype.TCellCD4Pos | hca7.cellType == Celltype.TCellCD8Pos | hca7.cellType == Celltype.TCellReg |  hca7.cellType == Celltype.TCell);
hca7_t.name = 'hca7_t';
hca7_b = hca7.cellSubset(hca7.cellType == Celltype.BCell);
hca7_b.name = 'hca7_b';
hca7_mon = hca7.cellSubset(hca7.cellType == Celltype.Monocyte);
hca7_mon.name = 'hca7_mon';
hca8_t = hca8.cellSubset(hca8.cellType == Celltype.TCellCD4Pos | hca8.cellType == Celltype.TCellCD8Pos | hca8.cellType == Celltype.TCellReg |  hca8.cellType == Celltype.TCell);
hca8_t.name = 'hca8_t';
hca8_b = hca8.cellSubset(hca8.cellType == Celltype.BCell);
hca8_b.name = 'hca8_b';
hca8_mon = hca8.cellSubset(hca8.cellType == Celltype.Monocyte);
hca8_mon.name = 'hca8_mon';

%% ovasc
ovasc = SCDep.scd_ovasc;
ovm = ovasc.cellSubset(ovasc.cellType == Celltype.MacrophageOrMonocyte);

%% PBMC68k
pbmc68000 = SCDep.scd_pbmc68000;
%{
figure
histogram(categorical(CelltypeId2CelltypeName(pbmc68000.cellType)));
title('PBMC68k cell types');
%}
t68000 = pbmc68000.cellSubset(pbmc68000.cellType == Celltype.TCellCD4Pos | pbmc68000.cellType == Celltype.TCellCD8Pos | pbmc68000.cellType == Celltype.TCellReg);
t68000.name = 't68000';
b68000 = pbmc68000.cellSubset(pbmc68000.cellType == Celltype.BCell);
b68000.name = 'b68000';
mon68000 = pbmc68000.cellSubset(pbmc68000.cellType == Celltype.Monocyte);
mon68000.name = 'mon68000';
nk68000 = pbmc68000.cellSubset(pbmc68000.cellType == Celltype.NKCell);
nk68000.name = 'nk68000';




%% Setup datasets and run DSAVE on all samples

%note: Removed the macrophages due to that there was only one sample in the
%breast cancer dataset that happened to be an outlier. Without it, the
%matrix did not have enough rank...
datasets = { ...
    lchp5t, ...
    lchp4t, ...
    lchp3t, ...
    lcp5t, ...
    lcp5mal, ...
    lcp5b, ...
    lcp3t, ...
    lcp3mal, ...
    bc1b_t, ...
    bc1n_t, ...
    bc1t_t, ...
    bc2ln_t, ...
    bc2ln_b, ...
    bc2t_t, ...
    bc4b_t, ...
    bc4b_b, ...
    bc4b_nk, ...
    bc4t_t, ...
    bc6t_t, ...
    bc8t_t, ...
    hca1_t, ...
    hca1_b, ...
    hca1_mon, ...
    hca1_nk, ...
    hca2_t, ...
    hca2_b, ...
    hca2_mon, ...
    hca2_nk, ...
    hca3_t, ...
    hca3_b, ...
    hca4_t, ...
    hca4_b, ...
    hca4_mon, ...
    hca4_nk, ...
    hca5_t, ...
    hca5_b, ...
    hca5_mon, ...
    hca5_nk, ...
    hca6_t, ...
    hca6_b, ...
    hca6_mon, ...
    hca7_t, ...
    hca7_b, ...
    hca7_mon, ...
    hca8_t, ...
    hca8_b, ...
    hca8_mon, ...
    t68000, ...
    b68000, ...
    mon68000, ...
    nk68000 ...
};

numSets = size(datasets,2);


%create a special template with only 1000 cells to allow for more
%populations
bc2t_ = SCDep.scd_bc2.cellSubset(SCDep.scd_bc2.cellType == Celltype.TCellCD4Pos | SCDep.scd_bc2.cellType == Celltype.TCellCD8Pos | SCDep.scd_bc2.cellType == Celltype.TCellReg);
bc2t_bc4tumor = bc2t_.cellSubset(strcmp(bc2t_.sampleIds, 'BC4_TUMOR'));
b10000 = SCDep.scd_pbmcb10000;
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

dss = {ovm,bc2t_bc4tumor, b10000, scd_GSE112845_cd8};
templSpec = DSAVEGenerateTemplateInfo(bc2t_bc4tumor, datasets, 1000, 750, 0.025, 0.025);

%check datasets
%{
for i = 1:numSets
    disp(strcat('dataset: ', num2str(i)));
    datasets{1,i}
    if (size(datasets{1,i}.cellIds,2) == 0)
        break;
    end
end
%}

scores = zeros(numSets,1);

for i = 1:numSets
    disp(strcat('running set: ', num2str(i)));
    DSAVERes = CalcDSAVE(datasets{1,i}, templSpec)
    scores(i,1) = DSAVERes.DSAVEScore;
end

%% setup design matrix and run regression
%intercept = t cell + hca + blood
%columns:
%icpt lc bc	oc	68k	mon	mal	b	mac	nk	tum	ht	ln
%{
dm = [ ...
1	1	0	0	0	0	0	0	0	0	0	1	0; ...	
1	1	0	0	0	0	0	0	1	0	0	1	0; ...	
1	1	0	0	0	0	0	0	0	0	0	1	0; ...	
1	1	0	0	0	0	0	0	0	0	0	1	0; ...	
1	1	0	0	0	0	0	0	0	0	1	0	0; ...	
1	1	0	0	0	0	1	0	0	0	1	0	0; ...	
1	1	0	0	0	0	0	1	0	0	1	0	0; ...	
1	1	0	0	0	0	0	0	0	0	1	0	0; ...	
1	1	0	0	0	0	1	0	0	0	1	0	0; ...	
1	1	0	0	0	0	0	0	1	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	1	0; ...	
1	0	1	0	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	0	1; ...	
1	0	1	0	0	0	0	1	0	0	0	0	1; ...	
1	0	1	0	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	0	0; ...	
1	0	1	0	0	0	0	1	0	0	0	0	0; ...	
1	0	1	0	0	0	0	0	0	1	0	0	0; ...	
1	0	1	0	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	1	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	1	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	1	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	1	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	1	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	1	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	1	0	0	0	0	1	0	1	0	0; ...	
1	0	0	0	1	0	0	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	1	0	0	0	0	0; ...	
1	0	0	0	1	1	0	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	0	0	1	0	0	0; ...	
];
%}

% intercept lc	bc	oc	68k	mon	mal	b	mac	nk	tml	tmb	htb	ln
%{
dm = [ ...
1	1	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	1	0	0	0	0	0	0	1	0	0	0	0	0; ...	
1	1	0	0	0	0	0	0	0	0	0	0	0   0; ...	
1	1	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	1	0	0	0	0	0	0	0	0	1	0	0	0; ...	
1	1	0	0	0	0	1	0	0	0	1	0	0	0; ...	
1	1	0	0	0	0	0	1	0	0	1	0	0	0; ...	
1	1	0	0	0	0	0	0	0	0	1	0	0	0; ...	
1	1	0	0	0	0	1	0	0	0	1	0	0	0; ...	
1	1	0	0	0	0	0	0	1	0	1	0	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	0	1	0; ...	
1	0	1	0	0	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	0	0	1; ...	
1	0	1	0	0	0	0	1	0	0	0	0	0	1; ...	
1	0	1	0	0	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	1	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	1	0	0	0	0	0	0	1	0	0	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	1	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	1	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	1	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	1	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	1	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	1	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0	0	0; ...	
1	0	0	1	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	1	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	1	1	0	0	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	0	0	1	0	0	0	0; ...	
];
%}

%icpt lc bc	68k	mon	b	nk	mal	tml	tmb	htb	ln
dm = [ ...
1	1	0	0	0	0	0	0	0	0	0	0; ...	
1	1	0	0	0	0	0	0	0	0	0   0; ...	
1	1	0	0	0	0	0	0	0	0	0	0; ...	
1	1	0	0	0	0	0	0	1	0	0	0; ...	
1	1	0	0	0	0	0	1	1	0	0	0; ...	
1	1	0	0	0	1	0	0	1	0	0	0; ...	
1	1	0	0	0	0	0	0	1	0	0	0; ...	
1	1	0	0	0	0	0	1	1	0	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	0; ...	
1	0	1	0	0	0	0	0	0	0	1	0; ...	
1	0	1	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	1; ...	
1	0	1	0	0	1	0	0	0	0	0	1; ...	
1	0	1	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	0	0	0; ...	
1	0	1	0	0	1	0	0	0	0	0	0; ...	
1	0	1	0	0	0	1	0	0	0	0	0; ...	
1	0	1	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	1	0	0; ...	
1	0	1	0	0	0	0	0	0	1	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	1	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	0	0	0	0	0	0	0	0	0; ...	
1	0	0	0	0	1	0	0	0	0	0	0; ...	
1	0	0	0	1	0	0	0	0	0	0	0; ...	
1	0	0	1	0	0	0	0	0	0	0	0; ...	
1	0	0	1	0	1	0	0	0	0	0	0; ...	
1	0	0	1	1	0	0	0	0	0	0	0; ...	
1	0	0	1	0	0	1	0	0	0	0	0; ...	
];


regrRes = dm\scores

%Export to text for analysis in R
headings = {'score', 'icpt', 'ds_lc', 'ds_bc', 'ds_68k', 'ct_mon', 'ct_b',	'ct_nk', 'ct_mal', 'ti_tumlung', 'ti_tumbreast', 'ti_healtbreast', 'ti_lymphnode'}
rowNames = cellfun(@(x) x.name, datasets, 'UniformOutput', false);
rowNames = rowNames.';
t = array2table([scores dm]);
t.Properties.RowNames = rowNames;
t.Properties.VariableNames = headings;
writetable(t,'../TempData/RegrData.txt','Delimiter','\t','WriteRowNames',true);


%intercept lc        bc         oc        68k       mon        mal       b        mac	    nk	      tum	    ht	       ln
%0.0271   -0.0567    0.1213    0.0256    0.0254    0.0015    0.0191   -0.0006   -0.0552    0.0004    0.0573    0.1022    0.1071

%int        lc        bc        oc       68k	   mon       mal        b        mac        nk       tml       tmb       htb	ln
%0.0290    0.0149    0.1208    0.0738    0.0253   -0.0004   -0.0070   -0.0045   -0.0479   -0.0014    0.0098    0.0368    0.2083  0.1076

%take a look at macrophages, they don't look right:
sel = dm(:,11) == 1;
scores(sel)
%test removing row 22, which is a low macrophage point
dm2 = dm;
scores2 = scores;
dm2(22,:) = [];
scores2(22,:) = [];
regrRes = dm2\scores2;
regrRes.'

%int        lc        bc        oc       68k	   mon       mal        b        mac        nk       tml       tmb       htb	ln
%0.0281    0.0038    0.1210    0.0268    0.0254    0.0005    0.0055   -0.0026    0.0000   -0.0005    0.0093    0.0517    0.2089 0.1074
