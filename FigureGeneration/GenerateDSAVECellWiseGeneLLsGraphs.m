    
%% Data Prep

%templInfo = DSAVEGetStandardTemplate();


[lc,lch] = SCDep.scd_lc;
lct = lc.cellSubset(lc.cellType == Celltype.TCellCD4Pos | lc.cellType == Celltype.TCellCD8Pos | lc.cellType == Celltype.TCellReg);
lctSub = lct.cellSubset(1:2000);%do not do random sampling for reproducability reasons.
lctSubTPM = TPM(lctSub); 
lcb = lc.cellSubset(lc.cellType == Celltype.BCell);
lcmast = lc.cellSubset(lc.cellType == Celltype.Mast);
lcpdc = lc.cellSubset(lc.cellType == Celltype.Dendritic);
lcer = lc.cellSubset(lc.cellType == Celltype.Erythroblast);
lchb = lch.cellSubset(lch.cellType == Celltype.BCell);
lcbSub = lcb.cellSubset(1:2000);%do not do random sampling for reproducability reasons.
lcbSubTPM = TPM(lcbSub); 


%% T cells

progbar = ProgrBar('Cell-wise: Fig E');

tic
[llsLct, genesLct, genellsLct] = DSAVEGetSingleCellDivergence(lctSub, 200, progbar.GetSubContext(1));
progbar.Done();
toc

%find the lowest ll value for each gene, and take the average of that:

%also need to remove only NaN cells if there are any - that is not the case
%here however
worstVals = nanmin(genellsLct, [], 2);
nanGenesSel = isnan(worstVals);
noNanGenes = genesLct(~nanGenesSel);
worstValsNoNan = worstVals(~nanGenesSel);
[~,worstValsIndices] = sort(worstValsNoNan);
worstValsNoNan(worstValsIndices(1:40))
noNanGenes(worstValsIndices(1:40))
%    {'HBB'      } - red blood cell precursor
%    {'IGKC'     } - Bcell/plasma cell
%    {'IGLC2'    } - Bcell/plasma cell
%    {'SLC52A2'  } - riboflavin (vitamin B2) transporter, cell 212, high brain and spinal cord cells
%    {'CPTP'     } - Unclear, cell 210, possibly brain tissue?
%    {'CCL4L2'   } - Cytokine, may be ok
%    {'IGHG1'    } - Bcell/plasma cell
%    {'EAPP'     } - Proliferation,cell 913, may be ok?
%    {'IGHA1'    } - Bcell/plasma cell
%    {'IGLC3'    } - Bcell/plasma cell
%    {'RPS15'    } - Ribosomal protein, unclear what role it plays here, cell 205 very high
%    {'C3AR1'    } - Exists in CD8+ T cells, may be ok
%    {'PITPNA-AS1'} - unclear
%    {'RSRP1'    } - unclear, a bit high in cell X
%    {'TMEM116'  } - unclear
%    {'RBMX'     } - unclear
%    {'STK17A'   } - Induces apoptosis, several cells with high expression, unclear what role this plays here
%    {'ATP6V1G1' } - ATPase related, unclear what role it plays here
%    {'HBA1'     } - red blood cell precursor
%    {'HBA2'     } - red blood cell precursor
%    {'S100A4'   } - unclear, exists in many cell types
%    {'DNAAF2'   } - unclear
%    {'PSMB8-AS1'} - non-coding, unclear
%    {'SLC39A1'  } - zink transporter, unclear
%    {'LDHA'      } - Lactate dehydrogenase A, often overexpressed in cancers
%    {'C16orf54' } - unclear
%    {'CCL4'     } - chemokine, produced by macrophages and CD8 T cells, ok
%    {'ANAPC16'   } - cell cycle-related, unclear
%    {'SAMSN1'   } - B cell suppressor, signature gene for TRegs (https://www.sciencedirect.com/science/article/pii/S1074761319300019), ok
%    {'DYNLT3'   } - unclear
%    {'MSANTD3'   } - unclear
%    {'RDH14'    } - unclear, metabolic gene
%    {'ZRANB3'   } - unclear
%    {'RAD21'     } - DNA repair, unclear
%    {'ORMDL2'    } - unclear
%    {'PPP1R7'   } - unclear
%    {'PPIB'     } - unclear
%    {'NFKBIA'    } - unclear
%    {'KLF6'      } - unclear
%    {'NAA16'     } - unclear



hbb = genellsLct(strcmp(genesLct,'HBB'),:)
hbb2 = lct.data(strcmp(lct.genes,'HBB'),:)

lctSubTPM.data(strcmp(lctSubTPM.genes,'LDHA'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'HBB'),:)

lctSubTPM.data(strcmp(lctSubTPM.genes,'IGLC2'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'IGLC3'),:)

lctSubTPM.data(strcmp(lctSubTPM.genes,'IGKC'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'SLC52A2'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'CPTP'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'EAPP'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'CCL4L2'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'NKG7'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'RPS15'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'STK17A'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'PTGDS'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'RSRP1'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'CCL4'),:)
lctSubTPM.data(strcmp(lctSubTPM.genes,'FTL'),:)

lcb.data(strcmp(lcb.genes,'IGKC'),:)

%plot against number of umis
numUMIsLct = sum(lctSub.data,1);

%find bad cells:
col = repmat([0.8,0.8,0.8], length(numUMIsLct),1);
threshold = 20000;

%red blood cell precursors:
selHBB = lctSubTPM.data(strcmp(lctSubTPM.genes,'HBB'),:) > threshold;
selHBA1 = lctSubTPM.data(strcmp(lctSubTPM.genes,'HBA1'),:) > threshold;
selHBA2 = lctSubTPM.data(strcmp(lctSubTPM.genes,'HBA2'),:) > threshold;

selRb = selHBB|selHBA1|selHBA2;
sum(selRb)
%b/plasma cells:
selIGKC = lctSubTPM.data(strcmp(lctSubTPM.genes,'IGKC'),:) > threshold;
selIGLC2 = lctSubTPM.data(strcmp(lctSubTPM.genes,'IGLC2'),:) > threshold;
selIGLC3 = lctSubTPM.data(strcmp(lctSubTPM.genes,'IGLC3'),:) > threshold;
selIGHA1 = lctSubTPM.data(strcmp(lctSubTPM.genes,'IGHA1'),:) > threshold;
selIGHG1 = lctSubTPM.data(strcmp(lctSubTPM.genes,'IGHG1'),:) > threshold;
selB = selIGKC | selIGLC2 | selIGLC3 | selIGHA1 | selIGHG1;

selLDHA = lctSubTPM.data(strcmp(lctSubTPM.genes,'LDHA'),:) > threshold;

falselyClassifiedT = selRb | selB;

group = zeros(2000,1);
group(selRb) = 1;
group(selB) = 2;
group(selLDHA) = 3;
groupT = group;
colT = [0.8 0.8 0.8; 1 0 0; 0 0.5 0; 0 0 0];

logNumUMIsLct = log2(numUMIsLct);

figure
gscatter(llsLct,logNumUMIsLct,group,colT,'o',3);
xlabel('Log likelihood')
ylabel('Log2(UMI counts)')
title('Classification of Divergent T Cells');
set(gca,'FontSize',11);
axis([-1000 -450 7 16]);

%make sure the colored plots are on top, so draw them again:
hold on;
scatter(llsLct(selLDHA),logNumUMIsLct(selLDHA),10,[0 0 0]);
hold on;
scatter(llsLct(selB),logNumUMIsLct(selB),10,[0,0.5,0]);
hold on;
scatter(llsLct(selRb),logNumUMIsLct(selRb),10,[1,0,0]);
legend({'No anomaly found', 'Red blood cell precursors', 'B/plasma cells', 'High lactate dehydrogenase A'});



%% Follicular B Cells

figure
histogram(lcb.subCellType)
folb = lcb.cellSubset(lcb.subCellType == 0 | lcb.subCellType == 1);
folSubTPM = TPM(folb);


progbar = ProgrBar('Cell-wise: Follicular B cells');

tic
[llsFol, genesFol, genellsFol] = DSAVEGetSingleCellDivergence(folb, 200, progbar.GetSubContext(1));
progbar.Done();
toc

%find the lowest ll value for each gene, and take the average of that:

%also need to remove only NaN cells if there are any - that is not the case
%here however
worstVals = nanmin(genellsFol, [], 2);
nanGenesSel = isnan(worstVals);
noNanGenes = genesFol(~nanGenesSel);
worstValsNoNan = worstVals(~nanGenesSel);
[~,worstValsIndices] = sort(worstValsNoNan);
worstValsNoNan(worstValsIndices(1:40))
noNanGenes(worstValsIndices(1:40))

%    {'IGHA1'   } - Immunoglobin
%    {'HBB'     } - red blood cell precursor
%    {'IGKC'    } - Immunoglobin
%    {'IGLC2'   } - Immunoglobin
%    {'IGLC3'   } - Immunoglobin
%    {'TOMM7'   }
%    {'MALAT1'  }
%    {'RPL35A'  }
%    {'COPE'    }
%    {'MED11'   }
%    {'SPP1'    }
%    {'IGHG3'   } - Immunoglobin
%    {'TACO1'   }
%    {'SCGB1A1' }
%    {'DOPEY2'  }
%    {'LTA4H'   }
%    {'MORC3'   }
%    {'SHMT2'   }
%    {'SYPL1'   }
%    {'MRPL34'  }
%    {'SCGB3A1' }
%    {'GNAS'    }
%    {'LAT2'    }
%    {'HSPA9'   }
%    {'PLP2'    }
%    {'RPL26'   }
%    {'SCP2'    }
%    {'CD27'    } - Mainly expressed on T cells
%    {'ERP29'   }
%    {'PIAS1'   }
%    {'IGHG1'   } - Immunoglobin
%    {'RPS2'    }
%    {'PHB'     }
%    {'SERGEF'  }
%    {'CLDND1'  }
%    {'C15orf40'}
%    {'HBA2'    } - red blood cell precursor
%    {'SLC25A6' }
%    {'KRCC1'   }
%    {'IGHG2'   } - Immunoglobin


folSubTPM.data(strcmp(folSubTPM.genes,'MALAT1'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'HBB'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'TPSB2'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'TPSAB1'),:)
lcmast.data(strcmp(lcmast.genes,'TPSB2'),:)
lcmast.data(strcmp(lcmast.genes,'TPSAB1'),:)
lcmast.data(strcmp(lcmast.genes,'MS4A2'),:)
lcmast.data(strcmp(lcmast.genes,'CD19'),:)
lcmast.data(strcmp(lcmast.genes,'IGHG1'),:)

lcpdc.data(strcmp(lcpdc.genes,'HBB'),:)

folSubTPM.data(strcmp(folSubTPM.genes,'SLC52A2'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'CPTP'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'EAPP'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'CCL4L2'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'NKG7'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'RPS15'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'STK17A'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'PTGDS'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'RSRP1'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'CCL4'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'FTL'),:)



%plot against number of umis
numUMIsFol = sum(folb.data,1);

%find bad cells:
threshold = 10000;

%red blood cell precursors:
selHBB = folSubTPM.data(strcmp(folSubTPM.genes,'HBB'),:) > threshold;
selHBA1 = folSubTPM.data(strcmp(folSubTPM.genes,'HBA1'),:) > threshold; %not among 40 most varying
selHBA2 = folSubTPM.data(strcmp(folSubTPM.genes,'HBA2'),:) > threshold; %not among 40 most varying

selRb = selHBB|selHBA1|selHBA2;
folSubTPM.data(strcmp(folSubTPM.genes,'HBB'),:)
sum(selRb)

%mast cells/basophils:
%    {'TPSB2' } - Serine protease, expressed in mast cells
%    {'TPSAB1'} - Serine protease, expressed in mast cells and basophils

selTPSB2 = folSubTPM.data(strcmp(folSubTPM.genes,'TPSB2'),:) > 5000;%the histogram looks separated here
selTPSB2Low = folSubTPM.data(strcmp(folSubTPM.genes,'TPSB2'),:) > 0;%presence of the gene
selTPSAB1 = folSubTPM.data(strcmp(folSubTPM.genes,'TPSAB1'),:) > 5000;%the histogram looks separated here
selBas = selTPSAB1 & ~selTPSB2;%TPSB2 is mainly expressed in mast cells, so if it has both, it is likely a mast cell
selMast = (selTPSB2 | selTPSAB1) & ~selBas;

%lung tissue:
selSFTPC = folSubTPM.data(strcmp(folSubTPM.genes,'SFTPC'),:) > 10000;
selLung = selSFTPC;
%SFTPCExpr = log2(folSubTPM.data(strcmp(folSubTPM.genes,'SFTPC'),:) + 1);
%figure
%histogram(SFTPCExpr)

TPSB2Expr = log2(folSubTPM.data(strcmp(folSubTPM.genes,'TPSB2'),:) + 1);
TPSAB1Expr = log2(folSubTPM.data(strcmp(folSubTPM.genes,'TPSAB1'),:) + 1);
figure
histogram(TPSB2Expr)

folSubTPM.data(strcmp(folSubTPM.genes,'IGHA1'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'IGKC'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'IGLC2'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'IGLC3'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'IGHG3'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'IGHG1'),:)
folSubTPM.data(strcmp(folSubTPM.genes,'IGHG2'),:)



selIGHA1 = folSubTPM.data(strcmp(folSubTPM.genes,'IGHA1'),:) > 15000;
selIGKC = folSubTPM.data(strcmp(folSubTPM.genes,'IGKC'),:) > 20000;%seems to be expressed a bit
selIGLC2 = folSubTPM.data(strcmp(folSubTPM.genes,'IGLC2'),:) > 15000;
selIGLC3 = folSubTPM.data(strcmp(folSubTPM.genes,'IGLC3'),:) > 15000;
selIGHG3 = folSubTPM.data(strcmp(folSubTPM.genes,'IGHG3'),:) > 20000;%seems to be expressed a bit
selIGHG1 = folSubTPM.data(strcmp(folSubTPM.genes,'IGHG1'),:) > 30000;%seems to be expressed a bit
selIGHG2 = folSubTPM.data(strcmp(folSubTPM.genes,'IGHG2'),:) > 15000;%seems to be expressed a bit

selIG = selIGHA1 | selIGKC | selIGLC2 | selIGLC3 | selIGHG3 | selIGHG1 | selIGHG2;
sum(selIG)

selSFTPC = folSubTPM.data(strcmp(folSubTPM.genes,'SFTPC'),:) > 10000;
selIG

%figure
%histogram(TPSAB1Expr)
%figure
%histogram(max(TPSB2Expr, TPSAB1Expr))

%selFTL = folSubTPM.data(strcmp(folSubTPM.genes,'FTL'),:) > threshold/2;

selT = folb.data(strcmp(folb.genes,'CD8A'),:) > 1;

falselyClassifiedFolb = selRb | selIG | selT;

selMan = [64 1088 845 864 1031 362];

group = zeros(size(folb.data,2),1);
group(selRb) = 1;
group(selIG) = 2;
group(selT) = 3;
%group(selMan) = 4;

%col = [0.8 0.8 0.8; 1 0 0; 0 0.5 0; 0 0 0; 0 0 1];
col = [0.8 0.8 0.8; 1 0 0; 0 0.5 0; 0 0 1];
logNumUMIsFol = log2(numUMIsFol);
figure
gscatter(llsFol,logNumUMIsFol,group,col,'o',3);
xlabel('Log likelihood')
ylabel('log2(UMI counts)')
title('Classification of Divergent Follicular B Cells');
set(gca,'FontSize',11);
axis([-2650 -800 8.6 14.5]);
%make sure the colored plots are on top, so draw them again:
hold on;
scatter(llsFol(selIG),logNumUMIsFol(selIG),10,[0,0.5,0]);
hold on;
scatter(llsFol(selRb),logNumUMIsFol(selRb),10,[1,0,0]);
hold on;
scatter(llsFol(selT),logNumUMIsFol(selT),10,[0,0,1]);
legend({'No anomaly found', 'Red blood cell precursors', 'Plasma cells', 'CD8+ T cells'});


%selMan = [64 1088 845 864 1031 362];

%look at highly divergent cells individually:
a = genellsFol(:,64);
[b,ind] = sort(a);
gi = genesFol(ind);
sel = ~isnan(b);
srt = [gi(sel) num2cell(b(sel))];
srt(1:50,:)

a = genellsFol(:,1088);
[b,ind] = sort(a);
gi = genesFol(ind);
sel = ~isnan(b);
srt = [gi(sel) num2cell(b(sel))];
srt(1:50,:)

%This one is likely a T cell or a doublet. Is variable in GZMB and
%expresses CD8A.
a = genellsFol(:,845);
[b,ind] = sort(a);
gi = genesFol(ind);
sel = ~isnan(b);
srt = [gi(sel) num2cell(b(sel))];
srt(1:50,:)
folSubTPM.data(strcmp(folSubTPM.genes,'CD8A'),:)
folb.data(strcmp(folSubTPM.genes,'CD8A'),:)


%Calculate DSAVE score before and after. Use a template of 1900 cells to
%make it possible to estimate the noise after removal of cells. Also don't
%remove any outliers; the outlier genes typically come from misclassified
%cells!

ovm = SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte);
bc2t = SCDep.scd_bc2.cellSubset(SCDep.scd_bc2.cellType == Celltype.TCellCD4Pos | SCDep.scd_bc2.cellType == Celltype.TCellCD8Pos | SCDep.scd_bc2.cellType == Celltype.TCellReg);
bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
b10000 = SCDep.scd_pbmcb10000;
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

datasets = {ovm,bc2t_bc4tumor, b10000, scd_GSE112845_cd8};
templ = DSAVEGenerateTemplateInfo(bc2t_bc4tumor, datasets, 1900, 750, 0, 0);

purT = lctSub.cellSubset(~falselyClassifiedT);
purfolb = folb.cellSubset(~falselyClassifiedFolb);

progbar = ProgrBar('Show improvement');
befScoreT = DSAVECalcBTMScore(lctSub, templ, progbar.GetSubContext(0.25));
aftScoreT = DSAVECalcBTMScore(purT, templ, progbar.GetSubContext(0.25));
befScoreB = DSAVECalcBTMScore(folb, templ, progbar.GetSubContext(0.25));
aftScoreB = DSAVECalcBTMScore(purfolb, templ, progbar.GetSubContext(0.25));

progbar.Done();

sprintf('TCells - before: %f, after: %f', befScoreT.DSAVEScore, aftScoreT.DSAVEScore)
sprintf('Fol B cells - before: %f, after: %f', befScoreB.DSAVEScore, aftScoreB.DSAVEScore)

%now, repeat this 100 times for each of those two populations, where we
%instead of removing the misclassified cells remove randomly selected cells
%Only need to run this code once, it takes ~4 hours, just use load below
randomT = zeros(100,1);
progbar = ProgrBar('Generate random T scores');
numCells = size(purT.data,2);
for i = 1:100
    score = DSAVECalcBTMScore(lctSub.randSample(numCells), templ, progbar.GetSubContext(0.01));
    randomT(i) = score.DSAVEScore;
end
figure
histogram(randomT);
save('../TempData/randomT.mat', 'randomT');



%Only need to run this code once, it takes ~4 hours, just use load below
randomFolB = zeros(100,1);
progbar = ProgrBar('Generate random FolB scores');
numCells = size(purfolb.data,2);
for i = 1:100
    score = DSAVECalcBTMScore(folb.randSample(numCells), templ, progbar.GetSubContext(0.01));
    randomFolB(i) = score.DSAVEScore;
end
figure
histogram(randomFolB);
save('../TempData/randomFolB.mat', 'randomFolB');

load('../TempData/randomT.mat');
[~,pT] = ttest(randomT, aftScoreT.DSAVEScore, 'Tail','right');
disp('P value for T cell DSAVE Score reduction: ');
pT % 4.1885e-86
load('../TempData/randomFolB.mat');
[~,pFolB] = ttest(randomFolB, aftScoreB.DSAVEScore, 'Tail','right');
disp('P value for Fol B cell DSAVE Score reduction: ');
pFolB % 3.7198e-91



%% Make a PCA plot

logSampT = LogTrans(lctSub.data,1);

[coeff, score, latent] = pca(logSampT.');
figure
gscatter(score(:,1),score(:,2),groupT,colT,'o',3);
xlabel('PC1')
ylabel('PC2')
title('PCA showing misclassified T Cells');
set(gca,'FontSize',11);
legend({'No anomaly found', 'Red blood cell precursors', 'B/plasma cells', 'High lactate dehydrogenase A'});
