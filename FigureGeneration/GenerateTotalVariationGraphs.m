%% Data initializatiom
totVarFromBulk = DSAVEGetPrecalcTotVarFromBulk()

%% Fig A

ub = 100000;
lb = 0.5;
n = 6000;

bulkMean1Vs1 = totVarFromBulk.totalVariation_1vs1(1,1);
bulkMean4Vs4 = totVarFromBulk.totalVariation_4vs4(1,1);

%create bulk data matrix
X = [0;repelem(n,99).'];
Y = repelem(bulkMean1Vs1,100).';
X = [X [0;repelem(n,99).']];
Y = [Y repelem(bulkMean4Vs4,100).'];


%set up datasets to process
[dsc,dsh] = SCDep.scd_lc();
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

dss = { SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte), ...
        SCDep.scd_livt, ...
        SCDep.scd_pbmc68000.cellSubset(SCDep.scd_pbmc68000.cellType == Celltype.TCellCD4Pos | SCDep.scd_pbmc68000.cellType == Celltype.TCellCD8Pos | SCDep.scd_pbmc68000.cellType == Celltype.TCellReg), ...
        dsc.cellSubset(dsc.cellType == Celltype.TCellCD4Pos | dsc.cellType == Celltype.TCellCD8Pos | dsc.cellType == Celltype.TCellReg), ...
        dsh.cellSubset(dsh.cellType == Celltype.TCellCD4Pos | dsh.cellType == Celltype.TCellCD8Pos | dsh.cellType == Celltype.TCellReg), ...
        scd_GSE112845_cd8 ...
      };
legends = { 'single bulk sample','mean of 4 bulk samples','OC - macrophages','LIVC - T cells', 'PBMC68k - T cells', 'LC - tumor T cells', 'LC - healthy tissue T cells', 'TCD8 - T cells'};
lineStyles = {'b--','k--','b-','k-','r-','m-','g-','c-'};

for i = 1:size(dss,2)
   vals = DSAVEGetTotalVariationVsPoolSize(dss{i}, n, ub, lb);
   X = [X vals(1,:).'];
   Y = [Y vals(2,:).'];
end


figure
%for some reason it doesn't work with the linestyles using a vector,
%so we'll have to loop
for i = 1:size(X,2)
    h = plot(X(:,i), Y(:,i), lineStyles{1,i});
    hold on
end
legend(legends);
xlabel('Pool size (number of cells)')
ylabel('Variation (R_m_e_a_n)')
title('Variation per Cell Pool Size. CPM > 0.5');
set(gca,'FontSize',11);





%% Fig B

ub = 2;
lb = 0.5;
n = 6000;

bulkMean1Vs1 = totVarFromBulk.totalVariation_1vs1(1,2);
bulkMean4Vs4 = totVarFromBulk.totalVariation_4vs4(1,2);

%create bulk data matrix
X2 = [0;repelem(n,99).'];
Y2 = repelem(bulkMean1Vs1,100).';
X2 = [X2 [0;repelem(n,99).']];
Y2 = [Y2 repelem(bulkMean4Vs4,100).'];


%set up datasets to process
[dsc,dsh] = SCDep.scd_lc();
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

dss = { SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte), ...
        SCDep.scd_livt, ...
        SCDep.scd_pbmc68000.cellSubset(SCDep.scd_pbmc68000.cellType == Celltype.TCellCD4Pos | SCDep.scd_pbmc68000.cellType == Celltype.TCellCD8Pos | SCDep.scd_pbmc68000.cellType == Celltype.TCellReg), ...
        dsc.cellSubset(dsc.cellType == Celltype.TCellCD4Pos | dsc.cellType == Celltype.TCellCD8Pos | dsc.cellType == Celltype.TCellReg), ...
        dsh.cellSubset(dsh.cellType == Celltype.TCellCD4Pos | dsh.cellType == Celltype.TCellCD8Pos | dsh.cellType == Celltype.TCellReg), ...
        scd_GSE112845_cd8 ...
      };

for i = 1:size(dss,2)
   vals = DSAVEGetTotalVariationVsPoolSize(dss{i}, n, ub, lb);
   X2 = [X2 vals(1,:).'];
   Y2 = [Y2 vals(2,:).'];
end


figure
%for some reason it doesn't work with the linestyles using a vector,
%so we'll have to loop
for i = 1:size(X2,2)
    h = plot(X2(:,i), Y2(:,i), lineStyles{1,i});
    hold on
end
legend(legends);
xlabel('Pool size (number of cells)')
ylabel('Variation (R_m_e_a_n)')
title('Variation per Cell Pool Size. 0.5 \leq CPM \leq 2');
set(gca,'FontSize',11);


%% Fig C


ub = 100000;
lb = 100;
n = 2000;

bulkMean1Vs1 = totVarFromBulk.totalVariation_1vs1(1,3);
bulkMean4Vs4 = totVarFromBulk.totalVariation_4vs4(1,3);

%create bulk data matrix
X3 = [0;repelem(n,99).'];
Y3 = repelem(bulkMean1Vs1,100).';
X3 = [X3 [0;repelem(n,99).']];
Y3 = [Y3 repelem(bulkMean4Vs4,100).'];


%set up datasets to process
[dsc,dsh] = SCDep.scd_lc();
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

dss = { SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte), ...
        SCDep.scd_livt, ...
        SCDep.scd_pbmc68000.cellSubset(SCDep.scd_pbmc68000.cellType == Celltype.TCellCD4Pos | SCDep.scd_pbmc68000.cellType == Celltype.TCellCD8Pos | SCDep.scd_pbmc68000.cellType == Celltype.TCellReg), ...
        dsc.cellSubset(dsc.cellType == Celltype.TCellCD4Pos | dsc.cellType == Celltype.TCellCD8Pos | dsc.cellType == Celltype.TCellReg), ...
        dsh.cellSubset(dsh.cellType == Celltype.TCellCD4Pos | dsh.cellType == Celltype.TCellCD8Pos | dsh.cellType == Celltype.TCellReg), ...
        scd_GSE112845_cd8 ...
      };

for i = 1:size(dss,2)
   vals = DSAVEGetTotalVariationVsPoolSize(dss{i}, n, ub, lb);
   X3 = [X3 vals(1,:).'];
   Y3 = [Y3 vals(2,:).'];
end


figure
%for some reason it doesn't work with the linestyles using a vector,
%so we'll have to loop
for i = 1:size(X3,2)
    h = plot(X3(:,i), Y3(:,i), lineStyles{1,i});
    hold on
end
legend(legends);
xlabel('Pool size (number of cells)')
ylabel('Variation (R_m_e_a_n)')
title('Variation per Cell Pool Size. CPM > 100');
set(gca,'FontSize',11);
axis([0 2000 0.038 0.52]);


%% Fig 4
%Typical histogram over gene expression

[dsc,dsh] = SCDep.scd_lc();
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

dss = { SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte), ...
        SCDep.scd_livt, ...
        SCDep.scd_pbmc68000.cellSubset(SCDep.scd_pbmc68000.cellType == Celltype.TCellCD4Pos | SCDep.scd_pbmc68000.cellType == Celltype.TCellCD8Pos | SCDep.scd_pbmc68000.cellType == Celltype.TCellReg), ...
        dsc.cellSubset(dsc.cellType == Celltype.TCellCD4Pos | dsc.cellType == Celltype.TCellCD8Pos | dsc.cellType == Celltype.TCellReg), ...
        dsh.cellSubset(dsh.cellType == Celltype.TCellCD4Pos | dsh.cellType == Celltype.TCellCD8Pos | dsh.cellType == Celltype.TCellReg), ...
        scd_GSE112845_cd8 ...
      };

legs = { 'OC macr.','LIVC T cells', 'PBMC68k T cells', 'LC tumor T cells', 'LC healthy tiss. T cells', 'TCD8 T cells'};
  
figure
for i = 1:size(dss, 2)
    avg = TPM(full(mean(dss{1,i}.data,2)));
    sel = avg >= 0.5 & avg <= 4000;

    d = log10(avg +0.05);
    a = linspace(-0.3,4,100);

    res = zeros(1, 100);
    for x = 1:size(a,2)
        s = d >= a(x)-0.1 & d <= a(x)+0.1;
        res(1,x) = sum(s);
    end
    res = res ./.2 ./ sum(sel);

    b = 10 .^a;

    semilogx(b.',res);
    hold on;
end

legend(legs);
xlabel('Gene expression (CPM)')
ylabel('Gene density')
ttl = 'Gene Density';
title(ttl);
axis([0.5 4000 0 0.8]);
set(gca,'FontSize',11);



