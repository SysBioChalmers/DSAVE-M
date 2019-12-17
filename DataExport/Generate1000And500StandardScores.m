%Generate templates

ovm = SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte);
bc2t = SCDep.scd_bc2.cellSubset(SCDep.scd_bc2.cellType == Celltype.TCellCD4Pos | SCDep.scd_bc2.cellType == Celltype.TCellCD8Pos | SCDep.scd_bc2.cellType == Celltype.TCellReg);
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
b10000 = SCDep.scd_pbmcb10000;
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

datasets = {ovm,bc2t_bc4tumor, b10000, scd_GSE112845_cd8};
templInfo1000 = DSAVEGenerateTemplateInfo(bc2t_bc4tumor, datasets, 1000, 750, 0.025, 0.025);
templInfo500 = DSAVEGenerateTemplateInfo(bc2t_bc4tumor, datasets, 500, 750, 0.025, 0.025);


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
resdata1000 = cell(1,numds);
scores1000 = zeros(1,numds);

progbar = ProgrBar('DSAVE Scores 1000:');


for i = 1:numds
    resdata1000{1,i} = DSAVECalcBTMScore(dss{1,i}, templInfo1000, progbar.GetSubContext(1/numds));
    scores1000(1,i) = resdata1000{1,i}.DSAVEScore;
end

progbar.Done();

disp('Scores 1000: Copy into R package');
scores1000 %just copy these values into the excel sheet

resdata500 = cell(1,numds);
scores500 = zeros(1,numds);

progbar = ProgrBar('DSAVE Scores 500:');


for i = 1:numds
    resdata500{1,i} = DSAVECalcBTMScore(dss{1,i}, templInfo500, progbar.GetSubContext(1/numds));
    scores500(1,i) = resdata500{1,i}.DSAVEScore;
end

progbar.Done();

disp('Scores 500: Copy into R package');
scores500 %just copy these values into the excel sheet

%export UMI distr and binning info

oldDir = SCDep.setPathToSource();
a = templInfo1000.UMIDistr;
save('../../ExportedData/templInfo1000UMIDistr.mat', 'a');
b = templInfo500.UMIDistr;
save('../../ExportedData/templInfo500UMIDistr.mat', 'b');
stdtempl = DSAVEGetStandardTemplate();
c = stdtempl.UMIDistr;
save('../../ExportedData/templInfo2000UMIDistr.mat', 'c');

a = templInfo1000.binningInfo;
save('../../ExportedData/templInfo1000BinningInfo.mat', 'a');
b = templInfo500.binningInfo;
save('../../ExportedData/templInfo500BinningInfo.mat', 'b');
c = stdtempl.binningInfo;
save('../../ExportedData/templInfo2000BinningInfo.mat', 'c');

SCDep.restoreDir(oldDir);


