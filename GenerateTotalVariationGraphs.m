%need to run SetupCTProfiles before running this file


%% Fig A

ub = 100000;
lb = 0.5;
n = 6000;

%calculate mean deviance between bulk samples (CD4+ T cells) from Blueprint
%read TMM normalization
samp = ImportTabSepSamples('../ImportableData/tcellCD4ProfilesTMMNormalized.txt');

bulkSamples = samp.cellSubset(32:39);
%remove all lowly expressed cells
means = mean(bulkSamples.data,2);
badGenes = means < lb | means > ub;
bulkSamples = bulkSamples.geneSubset(bulkSamples.genes(~badGenes)); 
numSamp = size(bulkSamples.sampleIds,2);
numGenes = size(bulkSamples.genes,1);
diffs = zeros(numGenes,numSamp*(numSamp-1)/2);
index = 1;

for i = 1:numSamp-1
   for j = i+1:numSamp
       diffs(:,index) = log2((bulkSamples.data(:,i)+0.05)./(bulkSamples.data(:,j)+0.05));
       index = index + 1;
   end
end

bulkMean1Vs1 = mean(mean(abs(diffs),2),1);

%a second bulk example, where the average of 4 samples is compared to the
%average of another 4.
%get all permutations of four
ind = 1:numSamp;
k = round(numSamp/2);
combs = nchoosek(ind,k);
numCombs = size(combs,1);
diffs = zeros(numGenes,numCombs);
for i = 1:numCombs
    a = bulkSamples.data(:,combs(i));
    notind = 1:numSamp;
    notind(combs(i,:)) = [];
    b = bulkSamples.data(:,notind);
    ma = mean(a,2);
    mb = mean(b,2);
    diffs(:,i) = log2((ma+0.05)./(mb+0.05));
end

bulkMean4Vs4 = mean(mean(abs(diffs),2),1);

%create bulk data matrix
X = [0;repelem(n,99).'];
Y = repelem(bulkMean1Vs1,100).';
X = [X [0;repelem(n,99).']];
Y = [Y repelem(bulkMean4Vs4,100).'];


%set up datasets to process
[dsc,dsh] = SCDep.scd_lc();
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

dss = { SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.paperClass == Celltype.MacrophageOrMonocyte), ...
        SCDep.scd_livt, ...
        SCDep.scd_pbmc68000.cellSubset(SCDep.scd_pbmc68000.paperClass == Celltype.TCellCD4Pos | SCDep.scd_pbmc68000.paperClass == Celltype.TCellCD8Pos | SCDep.scd_pbmc68000.paperClass == Celltype.TCellReg), ...
        dsc.cellSubset(dsc.paperClass == Celltype.TCellCD4Pos | dsc.paperClass == Celltype.TCellCD8Pos | dsc.paperClass == Celltype.TCellReg), ...
        dsh.cellSubset(dsh.paperClass == Celltype.TCellCD4Pos | dsh.paperClass == Celltype.TCellCD8Pos | dsh.paperClass == Celltype.TCellReg), ...
...%        SCDep.scd_pbmcb10000, ...
        scd_GSE112845_cd8 ...
      };
legends = { 'single bulk sample','mean of 4 bulk samples','ovarian cancer - macr.','liver cancer - T cells', 'blood - T cells', 'lung cancer - tumor T cells', 'lung cancer - healthy tissue T cells', 'blood cd8+ T cells'};
%lineStyles = {'--','-.','-.','-.','-.'};
lineStyles = {'b--','k--','b-','k-','r-','m-','g-','c-'};
%lineStyles = {'b--', 'b--', 'b--', 'b--'};
%lStyles = {'k','k--','k:','m','m--','m:','b','b--','b:','g','g--','g:','c','c--','c:','g-o','g--o','g:o','b-o','b--o','b:o','b-x','b--x','b:x','m-o','m--o','m:o','m-x','m--x','m:x'};

for i = 1:size(dss,2)
   vals = EvaluateClusterSize(dss{i}, n, ub, lb);
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
%h = plot(X, Y, 'linestyle',{'--','--','--','--'});
%hbc = get(h, 'Children');
%set(hbc{1}, 'FaceColor', 'r');
legend(legends);
xlabel('Pool size (number of cells)')
ylabel('Variation (R_m_e_a_n)')
title('Variation per Cell Pool Size. CPM > 0.5');
set(gca,'FontSize',11);





%% Fig B

ub = 2;
lb = 0.5;
n = 6000;

bulkSamples = samp.cellSubset(32:39);
%remove all lowly expressed cells
means = mean(bulkSamples.data,2);
badGenes = means < lb | means > ub;
bulkSamples = bulkSamples.geneSubset(bulkSamples.genes(~badGenes)); 
numSamp = size(bulkSamples.sampleIds,2);
numGenes = size(bulkSamples.genes,1);
diffs = zeros(numGenes,numSamp*(numSamp-1)/2);
index = 1;

for i = 1:numSamp-1
   for j = i+1:numSamp
       diffs(:,index) = log2((bulkSamples.data(:,i)+0.05)./(bulkSamples.data(:,j)+0.05));
       index = index + 1;
   end
end

bulkMean1Vs1 = mean(mean(abs(diffs),2),1);

%a second bulk example, where the average of 4 samples is compared to the
%average of another 4.
%get all permutations of four
ind = 1:numSamp;
k = round(numSamp/2);
combs = nchoosek(ind,k);
numCombs = size(combs,1);
diffs = zeros(numGenes,numCombs);
for i = 1:numCombs
    a = bulkSamples.data(:,combs(i));
    notind = 1:numSamp;
    notind(combs(i,:)) = [];
    b = bulkSamples.data(:,notind);
    ma = mean(a,2);
    mb = mean(b,2);
    diffs(:,i) = log2((ma+0.05)./(mb+0.05));
end

bulkMean4Vs4 = mean(mean(abs(diffs),2),1);


%create bulk data matrix
X2 = [0;repelem(n,99).'];
Y2 = repelem(bulkMean1Vs1,100).';
X2 = [X2 [0;repelem(n,99).']];
Y2 = [Y2 repelem(bulkMean4Vs4,100).'];


%set up datasets to process
[dsc,dsh] = SCDep.scd_lc();
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

dss = { SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.paperClass == Celltype.MacrophageOrMonocyte), ...
        SCDep.scd_livt, ...
        SCDep.scd_pbmc68000.cellSubset(SCDep.scd_pbmc68000.paperClass == Celltype.TCellCD4Pos | SCDep.scd_pbmc68000.paperClass == Celltype.TCellCD8Pos | SCDep.scd_pbmc68000.paperClass == Celltype.TCellReg), ...
        dsc.cellSubset(dsc.paperClass == Celltype.TCellCD4Pos | dsc.paperClass == Celltype.TCellCD8Pos | dsc.paperClass == Celltype.TCellReg), ...
        dsh.cellSubset(dsh.paperClass == Celltype.TCellCD4Pos | dsh.paperClass == Celltype.TCellCD8Pos | dsh.paperClass == Celltype.TCellReg), ...
...%        SCDep.scd_pbmcb10000, ...
        scd_GSE112845_cd8 ...
      };
%legends = { 'bulk 1:1','bulk 4:4','oc ascites - macrophages/monocytes','liver cancer - T cells', 'pbmc68000 - T cells', 'lung cancer - tumor T cells', 'lung cancer - healthy tissue T cells', 'GSE112845 - PBMC cd8+ T cells'};
%lineStyles = {'--','-.','-.','-.','-.'};
%lineStyles = {'b--','k--','b-','k-','r-','m-','g-','c-'};
%lineStyles = {'b--', 'b--', 'b--', 'b--'};
%lStyles = {'k','k--','k:','m','m--','m:','b','b--','b:','g','g--','g:','c','c--','c:','g-o','g--o','g:o','b-o','b--o','b:o','b-x','b--x','b:x','m-o','m--o','m:o','m-x','m--x','m:x'};

for i = 1:size(dss,2)
   vals = EvaluateClusterSize(dss{i}, n, ub, lb);
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
%h = plot(X, Y, 'linestyle',{'--','--','--','--'});
%hbc = get(h, 'Children');
%set(hbc{1}, 'FaceColor', 'r');
legend(legends);
xlabel('Pool size (number of cells)')
ylabel('Variation (R_m_e_a_n)')
title('Variation per Cell Pool Size. 0.5 \leq CPM \leq 2');
set(gca,'FontSize',11);


%% Fig C


ub = 100000;
lb = 100;
n = 2000;

bulkSamples = samp.cellSubset(32:39);
%remove all lowly expressed cells
means = mean(bulkSamples.data,2);
badGenes = means < lb | means > ub;
bulkSamples = bulkSamples.geneSubset(bulkSamples.genes(~badGenes)); 
numSamp = size(bulkSamples.sampleIds,2);
numGenes = size(bulkSamples.genes,1);
diffs = zeros(numGenes,numSamp*(numSamp-1)/2);
index = 1;

for i = 1:numSamp-1
   for j = i+1:numSamp
       diffs(:,index) = log2((bulkSamples.data(:,i)+0.05)./(bulkSamples.data(:,j)+0.05));
       index = index + 1;
   end
end

bulkMean1Vs1 = mean(mean(abs(diffs),2),1);

%a second bulk example, where the average of 4 samples is compared to the
%average of another 4.
%get all permutations of four
ind = 1:numSamp;
k = round(numSamp/2);
combs = nchoosek(ind,k);
numCombs = size(combs,1);
diffs = zeros(numGenes,numCombs);
for i = 1:numCombs
    a = bulkSamples.data(:,combs(i));
    notind = 1:numSamp;
    notind(combs(i,:)) = [];
    b = bulkSamples.data(:,notind);
    ma = mean(a,2);
    mb = mean(b,2);
    diffs(:,i) = log2((ma+0.05)./(mb+0.05));
end

bulkMean4Vs4 = mean(mean(abs(diffs),2),1);


%create bulk data matrix
X3 = [0;repelem(n,99).'];
Y3 = repelem(bulkMean1Vs1,100).';
X3 = [X3 [0;repelem(n,99).']];
Y3 = [Y3 repelem(bulkMean4Vs4,100).'];


%set up datasets to process
[dsc,dsh] = SCDep.scd_lc();
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

dss = { SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.paperClass == Celltype.MacrophageOrMonocyte), ...
        SCDep.scd_livt, ...
        SCDep.scd_pbmc68000.cellSubset(SCDep.scd_pbmc68000.paperClass == Celltype.TCellCD4Pos | SCDep.scd_pbmc68000.paperClass == Celltype.TCellCD8Pos | SCDep.scd_pbmc68000.paperClass == Celltype.TCellReg), ...
        dsc.cellSubset(dsc.paperClass == Celltype.TCellCD4Pos | dsc.paperClass == Celltype.TCellCD8Pos | dsc.paperClass == Celltype.TCellReg), ...
        dsh.cellSubset(dsh.paperClass == Celltype.TCellCD4Pos | dsh.paperClass == Celltype.TCellCD8Pos | dsh.paperClass == Celltype.TCellReg), ...
...%        SCDep.scd_pbmcb10000, ...
        scd_GSE112845_cd8 ...
      };
%legends = { 'bulk 1:1','bulk 4:4','oc ascites - macrophages/monocytes','liver cancer - T cells', 'pbmc68000 - T cells', 'lung cancer - tumor T cells', 'lung cancer - healthy tissue T cells', 'GSE112845 - PBMC cd8+ T cells'};
%lineStyles = {'--','-.','-.','-.','-.'};
%lineStyles = {'b--','k--','b-','k-','r-','m-','g-','c-'};
%lineStyles = {'b--', 'b--', 'b--', 'b--'};
%lStyles = {'k','k--','k:','m','m--','m:','b','b--','b:','g','g--','g:','c','c--','c:','g-o','g--o','g:o','b-o','b--o','b:o','b-x','b--x','b:x','m-o','m--o','m:o','m-x','m--x','m:x'};

for i = 1:size(dss,2)
   vals = EvaluateClusterSize(dss{i}, n, ub, lb);
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
%h = plot(X, Y, 'linestyle',{'--','--','--','--'});
%hbc = get(h, 'Children');
%set(hbc{1}, 'FaceColor', 'r');
legend(legends);
xlabel('Pool size (number of cells)')
ylabel('Variation (R_m_e_a_n)')
title('Variation per Cell Pool Size. CPM > 100');
set(gca,'FontSize',11);
axis([0 2000 0.038 0.52]);


%% Fig 4
%Typical histogram over gene expression

%some code for looking at the typical gene expression histogram

%a = linspace(1,3000,200);
%{
a = linspace(-0.3,3.2,100);

figure
for i = 1:size(dss,2)
    avg = TPM(full(mean(dss{1,i}.data,2)));
    d = log10(avg +0.05);
    res = zeros(1, 100);
    for x = 1:size(a,2)
        s = d >= a(1,x)-0.1 & d <= a(1,x)+0.1;
        res(1,x) = sum(s); 
        div = 10 .^ (a(1,x)+0.1) - 10 .^ (a(1,x)-0.1)
        res(1,x) = res(1,x) ./div;
        %res(1,x) = res(1,x) ./div ./ size(dss{1,i}.genes,1);
    end
    b = 10 .^a;
    %plot(b.',res);
    semilogy(b.',res);
    loglog(b.',res);
    hold on;
end

xlabel('Gene expression (CPM)')
ylabel('Genes per CPM')
title('Gene Expression Histograms');
legend(legends);
set(gca,'FontSize',11);
axis([0 100000 0.1 10000]);

%plot(a.',res);

%}


%Histogram over gene expression range
[dsc,dsh] = SCDep.scd_lc();
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

dss = { SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.paperClass == Celltype.MacrophageOrMonocyte), ...
        SCDep.scd_livt, ...
        SCDep.scd_pbmc68000.cellSubset(SCDep.scd_pbmc68000.paperClass == Celltype.TCellCD4Pos | SCDep.scd_pbmc68000.paperClass == Celltype.TCellCD8Pos | SCDep.scd_pbmc68000.paperClass == Celltype.TCellReg), ...
        dsc.cellSubset(dsc.paperClass == Celltype.TCellCD4Pos | dsc.paperClass == Celltype.TCellCD8Pos | dsc.paperClass == Celltype.TCellReg), ...
        dsh.cellSubset(dsh.paperClass == Celltype.TCellCD4Pos | dsh.paperClass == Celltype.TCellCD8Pos | dsh.paperClass == Celltype.TCellReg), ...
...%        SCDep.scd_pbmcb10000, ...
        scd_GSE112845_cd8 ...
      };

legs = { 'oc macr.','livc T cells', 'blood T cells', 'lc tumor T cells', 'lc healthy tiss. T cells', 'blood cd8+ T cells'};

  
figure
for i = 1:size(dss, 2);
    %some code for looking at the typical gene expression histogram
    avg = TPM(full(mean(dss{1,i}.data,2)));
    sel = avg >= 0.5 & avg <= 4000;
    %sum(sel)
    %figure
    %histogram(avg);
    %figure
    %histogram(log10(avg(sel)+.05));

    d = log10(avg +0.05);
    a = linspace(-0.3,4,100);

    res = zeros(1, 100);
    for x = 1:size(a,2)
        s = d >= a(x)-0.1 & d <= a(x)+0.1;
        res(1,x) = sum(s);%make this count as a 
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



