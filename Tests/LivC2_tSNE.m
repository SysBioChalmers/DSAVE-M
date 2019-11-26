%The TSNE plot here doesn't look random, the cell types partly cluster together. This shows that we get the
%cell classifications somewhat right.
ds = SCDep.scd_livc2().cellSubset(1:2000);

ds.data = full(ds.data);


A = tsne(ds.data.');

figure
gscatter(A(:,1),A(:,2),ds.cellType)
