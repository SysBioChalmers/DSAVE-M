%% Fig A - comparison between datasets

templInfo = DSAVEGetStandardTemplate();
hcacb = SCDep.scd_hca_cb;


%from breast cancer
bc2 = SCDep.scd_bc2;
bc2t = bc2.cellSubset(bc2.cellType == Celltype.TCellCD4Pos | bc2.cellType == Celltype.TCellCD8Pos | bc2.cellType == Celltype.TCellReg);
bc2t_bc1blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC1_BLOOD'));
bc2t_bc4blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_BLOOD'));
bc2t_2blood = bc2t_bc1blood.innerJoin(bc2t_bc4blood.randSample(size(bc2t_bc1blood.data,2)));%mix 50% from each patient

%ovarian cancer ascites
ovm = SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte);

%from lung cancer
[lc,lch] = SCDep.scd_lc;
lcht = lch.cellSubset(lch.cellType == Celltype.TCellCD4Pos | lch.cellType == Celltype.TCellCD8Pos | lch.cellType == Celltype.TCellReg);
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
t68000 = pbmc68000.cellSubset(pbmc68000.cellType == Celltype.TCellCD4Pos | pbmc68000.cellType == Celltype.TCellCD8Pos | pbmc68000.cellType == Celltype.TCellReg);

%b cells blood
b10000 = SCDep.scd_pbmcb10000;

%CD4+ memory t cells from blood
pbmctcd4mem10000 = SCDep.scd_pbmctcd4mem10000;

%from human cell atlas cord blood
hcat = hcacb.cellSubset(hcacb.cellType == Celltype.TCell | hcacb.cellType == Celltype.TCellCD4Pos | hcacb.cellType == Celltype.TCellCD8Pos);
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


numds = size(dss,2);
resdataA = cell(1,numds);
scoresA = zeros(1,numds);

progbar = ProgrBar('DSAVE Score 2: Fig A');


for i = 1:numds
    resdataA{1,i} = DSAVECalcBTMScore(dss{1,i}, templInfo, progbar.GetSubContext(1/numds));
    scoresA(1,i) = resdataA{1,i}.DSAVEScore;
end

progbar.Done();

disp('Comparison between datasets: Copy into excel sheet');
scoresA %just copy these values into the excel sheet

%% Fig 3E Suppl. (run fig A first)

progbar = ProgrBar('DSAVE Score 2: Fig 3E Suppl.');
resdata3ES = cell(1,numds);
scores3ES = zeros(1,numds);

for i = 1:numds
    resdata3ES{1,i} = DSAVECalcBTMScore(dss{1,i}, templInfo, progbar.GetSubContext(1/numds), false, 15, true);
    scores3ES(1,i) = resdata3ES{1,i}.DSAVEScore;
end

progbar.Done();

%figure out slope
b1 = (scores3ES.')\(scoresA.');

figure
plot([0 0.35], [0 b1*0.35]);
hold on
scatter(scores3ES, scoresA);
xlabel('BTM score using log transformation')
ylabel('BTM score')
title('BTM Score - Impact of Log Transformation');
%axis([0 1000 0 0.31]);
set(gca,'FontSize',11);

disp('Correlation matrix (3E Suppl.): ');
cc = corrcoef(scores3ES, scoresA)


%% Fig B - patient to patient variation

%from breast cancer
bc2 = SCDep.scd_bc2;
bc2t = bc2.cellSubset(bc2.cellType == Celltype.TCellCD4Pos | bc2.cellType == Celltype.TCellCD8Pos | bc2.cellType == Celltype.TCellReg);
bc2t_bc1tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC1_TUMOR'));
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
bc2t_2tumors = bc2t_bc1tumor.innerJoin(bc2t_bc4tumor.randSample(size(bc2t_bc1tumor.data,2)));%mix 50% from each patient
bc2t_bc1blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC1_BLOOD'));
bc2t_bc4blood = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_BLOOD'));
bc2t_2blood = bc2t_bc1blood.innerJoin(bc2t_bc4blood.randSample(size(bc2t_bc1blood.data,2)));%mix 50% from each patient

%from lung cancer
[lc,lch] = SCDep.scd_lc;
lct = lc.cellSubset(lc.cellType == Celltype.TCellCD4Pos | lc.cellType == Celltype.TCellCD8Pos | lc.cellType == Celltype.TCellReg);
lcht = lch.cellSubset(lch.cellType == Celltype.TCellCD4Pos | lch.cellType == Celltype.TCellCD8Pos | lch.cellType == Celltype.TCellReg);
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
hcat = hcacb.cellSubset(hcacb.cellType == Celltype.TCell | hcacb.cellType == Celltype.TCellCD4Pos | hcacb.cellType == Celltype.TCellCD8Pos);
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

%Generate a special template with only 1941 cells to be able to show both
%breast cancer patients
ovm = SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte);
bc2t = SCDep.scd_bc2.cellSubset(SCDep.scd_bc2.cellType == Celltype.TCellCD4Pos | SCDep.scd_bc2.cellType == Celltype.TCellCD8Pos | SCDep.scd_bc2.cellType == Celltype.TCellReg);
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
b10000 = SCDep.scd_pbmcb10000;
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

datasets = {ovm,bc2t_bc4tumor, b10000, scd_GSE112845_cd8};
templInfoSpec = DSAVEGenerateTemplateInfo(bc2t_bc4tumor, datasets, 1941, 570, 0.025, 0.025);

progbar = ProgrBar('DSAVE Score 2: Fig B');

for i = 1:numds
    resdata{1,i} = DSAVECalcBTMScore(dss{1,i}, templInfoSpec, progbar.GetSubContext(1/numds));
    scores(1,i) = resdata{1,i}.DSAVEScore;
end

progbar.Done();

disp('Patient to patient: Copy into excel sheet');
scores %just copy these values into the excel sheet

%% Fig C - effect of mixing cell types

%create special template with 1346 cells
ovm = SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte);
bc2t = SCDep.scd_bc2.cellSubset(SCDep.scd_bc2.cellType == Celltype.TCellCD4Pos | SCDep.scd_bc2.cellType == Celltype.TCellCD8Pos | SCDep.scd_bc2.cellType == Celltype.TCellReg);
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
b10000 = SCDep.scd_pbmcb10000;
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;
datasets = {ovm,bc2t_bc4tumor, b10000, scd_GSE112845_cd8};
spec1346Templ = DSAVEGenerateTemplateInfo(bc2t_bc4tumor, datasets, 1346, 750, 0.025, 0.025);




%from breast cancer
bc2 = SCDep.scd_bc2;
bc2_bc4blood = bc2.cellSubset(strcmp(bc2.sampleIds, 'BC4_BLOOD'));
bc2_bc4blood_t = bc2_bc4blood.cellSubset(bc2_bc4blood.cellType == Celltype.TCellCD4Pos | bc2_bc4blood.cellType == Celltype.TCellCD8Pos | bc2_bc4blood.cellType == Celltype.TCellReg);
bc2_bc4blood_b = bc2_bc4blood.cellSubset(bc2_bc4blood.cellType == Celltype.BCell);
numb = size(bc2_bc4blood_b.data,2);
bc2_bc4blood_mixbc50_50 = bc2_bc4blood_b.innerJoin(bc2_bc4blood_t.randSample(numb));%mix 50% from each patient

%from human cell atlas cord blood
hcacb = SCDep.scd_hca_cb;
hcat = hcacb.cellSubset(hcacb.cellType == Celltype.TCell | hcacb.cellType == Celltype.TCellCD4Pos | hcacb.cellType == Celltype.TCellCD8Pos);
hcab = hcacb.cellSubset(hcacb.cellType == Celltype.BCell);
hcat_cb1 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB1'));
hcab_cb1 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB1'));
num = size(hcat_cb1.data,2);
numb1 = size(hcab_cb1.data,2);
hcat_1bt50_50 = hcab_cb1.innerJoin(hcat_cb1.randSample(numb1));%50/50 b and t

%from lc
[lc,~] = SCDep.scd_lc;
lcpat5 = lc.cellSubset(strcmp(lc.sampleIds,'5'));
lcp5t = lcpat5.cellSubset(lcpat5.cellType == Celltype.TCellCD4Pos | lcpat5.cellType == Celltype.TCellCD8Pos | lcpat5.cellType == Celltype.TCellReg |  lcpat5.cellType == Celltype.TCell);
lcp5b = lcpat5.cellSubset(lcpat5.cellType == Celltype.BCell);
numb = size(lcp5b.data,2);
lc50_50 = lcp5b.innerJoin(lcp5t.randSample(numb));%50/50 b and t


dss = { bc2_bc4blood_t, ...
        bc2_bc4blood_b, ...
        bc2_bc4blood_mixbc50_50, ...
        lcp5t, ...
        lcp5b, ...
        lc50_50 ...
        hcat_cb1, ...
        hcab_cb1, ...
        hcat_1bt50_50 ...
      };

numds = size(dss,2);
resdata = cell(1,numds);
scores = zeros(1,numds);

progbar = ProgrBar('DSAVE Score 2: Fig C');

for i = 1:numds
    resdata{1,i} = DSAVECalcBTMScore(dss{1,i}, spec1346Templ, progbar.GetSubContext(1/numds));
    scores(1,i) = resdata{1,i}.DSAVEScore;
end

progbar.Done();

disp('Mixing cell types: Copy into excel sheet');
scores %just copy these values into the excel sheet


%% Supplementary 3A - technical validation using generated cell populations with added noise

templInfo = DSAVEGetStandardTemplate();

progbar = ProgrBar('DSAVE Score 2: Fig 3A Suppl.');

ovasc = SCDep.scd_ovasc;
b10000 = SCDep.scd_pbmcb10000;
pbmc68000 = SCDep.scd_pbmc68000;
bc2 = SCDep.scd_bc2;
hcacb = SCDep.scd_hca_cb;

bc2t = bc2.cellSubset(bc2.cellType == Celltype.TCellCD4Pos | bc2.cellType == Celltype.TCellCD8Pos | bc2.cellType == Celltype.TCellReg);
hcat = hcacb.cellSubset(hcacb.cellType == Celltype.TCell | hcacb.cellType == Celltype.TCellCD4Pos | hcacb.cellType == Celltype.TCellCD8Pos);
ovm = ovasc.cellSubset(ovasc.cellType == Celltype.MacrophageOrMonocyte);
t68000 = pbmc68000.cellSubset(pbmc68000.cellType == Celltype.TCellCD4Pos | pbmc68000.cellType == Celltype.TCellCD8Pos | pbmc68000.cellType == Celltype.TCellReg);

%throw away some cells to save computational time
hcatc = hcat.randSample(10000);
bc2tc = bc2t.randSample(10000);
t68000c = t68000.randSample(10000);

bc2tGenNoise = DSAVEGenerateSNODataset(bc2tc, progbar.GetSubContext(0.01), size(bc2tc.data, 2),0.7);
hcatGenNoise = DSAVEGenerateSNODataset(hcatc, progbar.GetSubContext(0.01), size(hcatc.data, 2),0.9);
ovmGenNoise = DSAVEGenerateSNODataset(ovm, progbar.GetSubContext(0.01), size(ovm.data, 2),1);
t68000GenNoise = DSAVEGenerateSNODataset(t68000c, progbar.GetSubContext(0.01), size(t68000c.data, 2),1.5);
b10000GenNoise = DSAVEGenerateSNODataset(b10000, progbar.GetSubContext(0.01), size(b10000.data, 2),2);

bc2tGenNoise1 = DSAVEGenerateSNODataset(bc2tc, progbar.GetSubContext(0.01), size(bc2tc.data, 2),1);
hcatGenNoise1 = DSAVEGenerateSNODataset(hcatc, progbar.GetSubContext(0.01), size(hcatc.data, 2),1);
t68000GenNoise1 = DSAVEGenerateSNODataset(t68000c, progbar.GetSubContext(0.01), size(t68000c.data, 2),1);
b10000GenNoise1 = DSAVEGenerateSNODataset(b10000, progbar.GetSubContext(0.01), size(b10000.data, 2),1);



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

numds = size(dss,2);
resdata = cell(1,numds);
scores = zeros(1,numds);

for i = 1:numds
    resdata{1,i} = DSAVECalcBTMScore(dss{1,i}, templInfo, progbar.GetSubContext(1/numds*0.9));
    scores(1,i) = resdata{1,i}.DSAVEScore;
end

progbar.Done();

disp('Technical validation: Copy into excel sheet');
scores %just copy these values into the excel sheet

%% Fig 3B Suppl. - Mix of monocytes and T cells
hcacb = SCDep.scd_hca_cb;
hca_cb1 = hcacb.cellSubset(strcmp(hcacb.sampleIds,'CB1'));
hcat = hca_cb1.cellSubset(hca_cb1.cellType == Celltype.TCell | hca_cb1.cellType == Celltype.TCellCD4Pos | hca_cb1.cellType == Celltype.TCellCD8Pos);
hcam = hca_cb1.cellSubset(hca_cb1.cellType == Celltype.Monocyte);
hcat2500 = hcat.randSample(2500);
hcam2500 = hcam.randSample(2500);
templInfo = DSAVEGetStandardTemplate();


progbar = ProgrBar('DSAVE Score 2: Fig 3B Suppl.');

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
    score = DSAVECalcBTMScore(ds, templInfo, progbar.GetSubContext(1/21));
    vals(1,i) = score.DSAVEScore;
end

progbar.Done();

figure
plot(xes,vals);
xlabel('Percent monocytes in mix');
ylabel('DSAVE variation score');
title('Variation per Cell Type Fractions');
set(gca,'FontSize',11);
axis([0 100 0 0.06]);

%% Fig 3D Suppl. - Evaluate number of iterations for the DSAVE score

ds = bc2t_2blood;

runs = 15;
iterations = 30;

resdata = cell(1,runs);
differenceCVs = cell(1,runs);

templInfo = DSAVEGetStandardTemplate();

progbar = ProgrBar('DSAVE Score 2: Fig 3D Suppl.');

for r = 1:runs
    [resdata{1,r}, differenceCVs{1,r}] = DSAVECalcBTMScore(ds, templInfo, progbar.GetSubContext(0.49/runs), false, iterations);
end

xvals = 1:iterations;
allScores = zeros(runs, iterations);
stds = zeros(1,iterations);
for r = 1:runs
    diffCVsRun = differenceCVs{1,r};
    allScores(r,:) = mean(diffCVsRun, 2);% This is not exactly the way it is done in the score, but since we're only using means, that should be fine
end



% now loop through all iteration sizes s and use the 1:s first iterations
% to form a mean. Then look at the std for those
for i = xvals
    meanVals = mean(allScores(:,1:i),2);
    stds(1,i) = std(meanVals, [], 1);
end

progbar.Done();


figure
plot(xvals, stds);
xlabel('Iterations');
ylabel('Standard deviation of DSAVE score');
title('Number of Iterations for the DSAVE Score');
set(gca,'FontSize',11);
axis([0 30 0 0.003]);

%repeat for a different template with fewer UMIs per cell
%create a special template with only 1000 cells to allow for more
%populations
ovm = SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte);
bc2t_ = SCDep.scd_bc2.cellSubset(SCDep.scd_bc2.cellType == Celltype.TCellCD4Pos | SCDep.scd_bc2.cellType == Celltype.TCellCD8Pos | SCDep.scd_bc2.cellType == Celltype.TCellReg);
bc2t_bc4tumor = bc2t_.cellSubset(strcmp(bc2t_.sampleIds, 'BC4_TUMOR'));
b10000 = SCDep.scd_pbmcb10000;
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

dss = {ovm,bc2t_bc4tumor, b10000, scd_GSE112845_cd8};
templSpec = DSAVEGenerateTemplateInfo(bc2t_bc4tumor, datasets, 1000, 750, 0.025, 0.025, progbar.GetSubContext(0.02));


resdata2 = cell(1,runs);
differenceCVs2 = cell(1,runs);

for r = 1:runs
    [resdata2{1,r}, differenceCVs2{1,r}] = DSAVECalcBTMScore(ds, templSpec, progbar.GetSubContext(0.49/runs), false, iterations);
end

allScores = zeros(runs, iterations);
stds2 = zeros(1,iterations);
for r = 1:runs
    diffCVsRun = differenceCVs2{1,r};
    allScores(r,:) = mean(diffCVsRun, 2);% This is not exactly the way it is done in the score, but since we're only using means, that should be fine
end

% now loop through all iteration sizes s and use the 1:s first iterations
% to form a mean. Then look at the std for those
for i = xvals
    meanVals = mean(allScores(:,1:i),2);
    stds2(1,i) = std(meanVals, [], 1);
end

hold on
plot(xvals, stds2);
legend({'Standard template - 2,000 cells', 'Specialized template - 1,000 cells'});
 
%% Test to create a PCA plot

[lct, lch] = SCDep.scd_lc();
[~,~,gse112845cd8] = SCDep.scd_GSE112845();
dsList = {SCDep.scd_hca_cb, lct, lch, SCDep.scd_pbmc68000(), SCDep.scd_pbmctcd4mem10000(), SCDep.scd_pbmcb10000(), gse112845cd8, SCDep.scd_bc2()};

dsList = SynchronizeGenes(dsList, [], true);

hca = dsList{1};
lct = dsList{2};
lch = dsList{3};
pbmc68k = dsList{4};
tcd4mem = dsList{5};
b10k = dsList{6};
gse112845cd8 = dsList{7};
bc = dsList{8};

%% HCA TCells and BCells from 8 patients, in total 16 samples
hcat = hca.cellSubset(hca.cellType == Celltype.TCellCD4Pos | hca.cellType == Celltype.TCellCD8Pos | hca.cellType == Celltype.TCellReg |  hca.cellType == Celltype.TCell );
hcab = hca.cellSubset(hca.cellType == Celltype.BCell);

hcat1 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB1'));
hcat2 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB2'));
hcat3 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB3'));
hcat4 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB4'));
hcat5 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB5'));
hcat6 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB6'));
hcat7 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB7'));
hcat8 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB8'));

hcab1 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB1'));
hcab2 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB2'));
hcab3 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB3'));
hcab4 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB4'));
hcab5 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB5'));
hcab6 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB6'));
hcab7 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB7'));
hcab8 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB8'));

%% PBMC68k : one sample of each, this is just one patient
pbmc68kt = pbmc68k.cellSubset(pbmc68k.cellType == Celltype.TCellCD4Pos | pbmc68k.cellType == Celltype.TCellCD8Pos | pbmc68k.cellType == Celltype.TCellReg |  pbmc68k.cellType == Celltype.TCell );
pbmc68kb = pbmc68k.cellSubset(pbmc68k.cellType == Celltype.BCell);

%% BC
bc1 = bc.cellSubset(strcmp(bc.sampleIds, 'BC1_BLOOD'));
bc4 = bc.cellSubset(strcmp(bc.sampleIds, 'BC4_BLOOD'));
bc1b = bc1.cellSubset(bc1.cellType == Celltype.BCell);%too few cells
bc1t = bc1.cellSubset(bc1.cellType == Celltype.TCellCD4Pos | bc1.cellType == Celltype.TCellCD8Pos | bc1.cellType == Celltype.TCellReg |  bc1.cellType == Celltype.TCell);
bc4b = bc4.cellSubset(bc4.cellType == Celltype.BCell);%1300 cells, a bit too few
bc4t = bc4.cellSubset(bc4.cellType == Celltype.TCellCD4Pos | bc4.cellType == Celltype.TCellCD8Pos | bc4.cellType == Celltype.TCellReg |  bc4.cellType == Celltype.TCell);

unique(bc.sampleIds)
figure
histogram(categorical(bc.sampleIds))


%% LC : use pat 3, 4 and 5 from the cancer and all patients for the healthy
%tissue; there are fewer cells there

lctpat3 = lct.cellSubset(strcmp(lct.sampleIds,'3'));
lctpat4 = lct.cellSubset(strcmp(lct.sampleIds,'4'));
lctpat5 = lct.cellSubset(strcmp(lct.sampleIds,'5'));

lctpat3t = lctpat3.cellSubset(lctpat3.cellType == Celltype.TCellCD4Pos | lctpat3.cellType == Celltype.TCellCD8Pos | lctpat3.cellType == Celltype.TCellReg |  lctpat3.cellType == Celltype.TCell );
lctpat4t = lctpat4.cellSubset(lctpat4.cellType == Celltype.TCellCD4Pos | lctpat4.cellType == Celltype.TCellCD8Pos | lctpat4.cellType == Celltype.TCellReg |  lctpat4.cellType == Celltype.TCell );
lctpat5t = lctpat5.cellSubset(lctpat5.cellType == Celltype.TCellCD4Pos | lctpat5.cellType == Celltype.TCellCD8Pos | lctpat5.cellType == Celltype.TCellReg |  lctpat5.cellType == Celltype.TCell );

lctpat3b = lctpat3.cellSubset(lctpat3.cellType == Celltype.BCell);
lctpat4b = lctpat4.cellSubset(lctpat4.cellType == Celltype.BCell);
lctpat5b = lctpat5.cellSubset(lctpat5.cellType == Celltype.BCell);

lcht = lch.cellSubset(lch.cellType == Celltype.TCellCD4Pos | lch.cellType == Celltype.TCellCD8Pos | lch.cellType == Celltype.TCellReg |  lch.cellType == Celltype.TCell );
lchb = lch.cellSubset(lch.cellType == Celltype.BCell);


%% TCD4Mem


%% B10k


%% GSE112845

%% Export them all to a Samples object
%dss = { hcat1, hcat2, hcat3, hcat4, hcat5, hcat6, hcat7, hcat8, hcab1, hcab2, hcab3, hcab4, hcab5, hcab6, hcab7, hcab8, ...
%        lctpat3t, lctpat4t, lctpat5t, lctpat3b, lctpat4b, lctpat5b, lcht, lchb, pbmc68kt, pbmc68kb, tcd4mem, b10k, gse112845cd8, };
dss = { hcat1, hcat2, hcat3, hcat4, hcat5, hcat6, hcat7, hcat8, ...
        lcht, pbmc68kt, tcd4mem, gse112845cd8, bc1t, bc4t };
numSets = size(dss,2);
numGenes = size(hcat1.data,1);

samp = Samples;
samp.data = zeros(numGenes,numSets);
%samp.sampleIds = {'hcat1', 'hcat2', 'hcat3', 'hcat4', 'hcat5', 'hcat6', 'hcat7', 'hcat8', 'hcab1', 'hcab2', 'hcab3', 'hcab4', 'hcab5', ...
%                  'hcab6', 'hcab7', 'hcab8', 'lctpat3t', 'lctpat4t', 'lctpat5t', 'lctpat3b', 'lctpat4b', 'lctpat5b', 'lcht', 'lchb', ...
%                  'pbmc68kt', 'pbmc68kb', 'tcd4mem', 'b10k', 'gse112845cd8', };
samp.sampleIds = {'hcat1', 'hcat2', 'hcat3', 'hcat4', 'hcat5', 'hcat6', 'hcat7', 'hcat8', ...
                  'lcht', 'pbmc68kt', 'tcd4mem', 'gse112845cd8', 'bc1t', 'bc4t'};
samp.genes = hcat1.genes;


for i = 1:numSets
    ds = dss{1,i};
    %we do not tpm normalize until after summing up all cells, since we
    %want the counts to truly represent the noise. This should not matter
    %much
    datasum = sum(ds.data,2);
    samp.data(:,i) = TPM(datasum);
end


%log transform
logSamp = LogTrans(samp.data,1);

[coeff, score, latent] = pca(logSamp.');
colors = [0 0 0 0 0 0 0 0 1 2 3 4 5 5];
figure
gscatter(score(:,1),score(:,2),colors);


%dlmwrite('scTotCounts.txt',totcounts,'\t');


