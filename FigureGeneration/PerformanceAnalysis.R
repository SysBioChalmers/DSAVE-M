#This file assumes that the DSAVE package is loaded (using open the DSAVE-R package, run devtools::load_all("."))
#We measure the execution time of the DSAVE R package - generates Fig A in S1 Note

library ("tictoc")
library ("Seurat")
#devtools::load_all(".")
extrDir <- downloadData("https://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz", "B10k")
dataDir = paste0(extrDir,"/filtered_matrices_mex/hg19")
bcells <- Read10X(data.dir = dataDir)

b1k = bcells[, sample(dim(bcells)[2], 1000)]
b5k = bcells[, sample(dim(bcells)[2], 5000)]
b10k = bcells[, sample(dim(bcells)[2], 10000)]

templ500 = DSAVEGetStandardTemplate500()
templ1000 = DSAVEGetStandardTemplate1000()
templ2000 = DSAVEGetStandardTemplate()

tb1k500s = rep(0,10)
tb5k500s = rep(0,10)
tb10k500s = rep(0,10)
tb10k1000s = rep(0,10)
tb10k2000s = rep(0,10)


for (i in 1:10) {
  tic()
  tmp = DSAVECalcBTMScore(b1k, templ500)
  t = toc()
  tb1k500s[i] = t$toc - t$tic

  tic()
  tmp = DSAVECalcBTMScore(b5k, templ500)
  t = toc()
  tb5k500s[i] = t$toc - t$tic

  tic()
  tmp = DSAVECalcBTMScore(b10k, templ500)
  t = toc()
  tb10k500s[i] = t$toc - t$tic

  tic()
  tmp = DSAVECalcBTMScore(b10k, templ1000)
  t = toc()
  tb10k1000s[i] = t$toc - t$tic

  tic()
  tmp = DSAVECalcBTMScore(b10k, templ2000)
  t = toc()
  tb10k2000s[i] = t$toc - t$tic
}

fnPerfScore="C:/Work/MatlabCode/components/DSAVE-M/TempData/perfScore.RData"
save(tb1k500s, tb5k500s, tb10k500s, tb10k1000s, tb10k2000s, file=fnPerfScore)
load(file=fnPerfScore)

tb1k500m = mean(tb1k500s)
tb5k500m = mean(tb5k500s)
tb10k500m = mean(tb10k500s)
tb10k1000m = mean(tb10k1000s)
tb10k2000m = mean(tb10k2000s)


#plot number of cells vs time
scoreX = c(1000, 5000, 10000)
scoreY = c(tb1k500m, tb5k500m, tb10k500m)
dfScore = data.frame(x=scoreX, y=scoreY)

library(ggplot2)

p1 = ggplot(data=dfScore, aes(x=x, y=y)) +
  #geom_bar(stat="identity", width = 1000) +
  geom_line() +
  geom_point() +
  labs(title="BTM Score: Cells",
       x ="Number of cells", y = "Time(s)") +
  scale_x_continuous(name ="Number of cells", breaks=c(1000,5000,10000)) +
  ylim(0,120) +
  theme_bw()
p1

#plot number of cells in template vs time
scoreX = c(500, 1000, 2000)
scoreY = c(tb10k500m, tb10k1000m, tb10k2000m)
dfScore = data.frame(x=scoreX, y=scoreY)

p2 = ggplot(data=dfScore, aes(x=x, y=y)) +
  geom_line() +
  geom_point() +
  #geom_bar(stat="identity", width = 400) +
  labs(title="BTM Score: Cells in template",
       x ="Number of cells in template", y = "Time(s)")  +
  scale_x_continuous(name ="Number of cells in template", breaks=c(500,1000,2000)) +
  ylim(0,600) +
  theme_bw()
p2


#also test the performance on pbmc68k, the whole population
extrDir <- downloadData("https://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz", "PBMC68k")
dataDir = paste0(extrDir,"/filtered_matrices_mex/hg19")
pbmc <- Read10X(data.dir = dataDir)

tic()
tmp = DSAVECalcBTMScore(pbmc, templ2000)
t = toc() #419.05, CPU at ~33%


#Cell wise vatiation:

cwNumCells = c(500, 1000, 2000, 4000, 10000)
numPoints = length(cwNumCells)
cw = matrix(0,nrow=numPoints, ncol=10)

for (iter in 1:10) {
  for (i in 1:numPoints) {
    sub = bcells[, sample(dim(bcells)[2], cwNumCells[i])]
    tic()
    div = DSAVEGetSingleCellDivergence(sub, minUMIsPerCell = 200)
    t = toc()
    cw[i,iter] = t$toc - t$tic
    gc()
  }
}

#runs at ~40% CPU

fnPerfCw="C:/Work/MatlabCode/components/DSAVE-M/TempData/perfCw.RData"
save(cw, file=fnPerfCw)
load(file=fnPerfCw)

cwVals = rowMeans(cw)

#plot number of cells in template vs time

dfCw = data.frame(x=cwNumCells, y=cwVals)

p3 = ggplot(data=dfCw, aes(x=x, y=y)) +
  geom_line() +
  geom_point() +
  #geom_bar(stat="identity", width = 400) +
  labs(title="Cell Divergence: Cells",
       x ="Number of cells", y = "Time(s)")  +
  #scale_x_continuous(name ="Number of cells", breaks=cwNumCells) +
  ylim(0,2500) +
  theme_bw()
p3

library(ggpubr)

figA = ggarrange(p1, p2, p3, nrow=2, ncol=2,labels=c("I","II","III"))
figA

fig__path = "Z:/projects/Partitioning of Cell-to-Cell Variation in Single-Cell RNA-Seq Data/paper/PLOS ONE/Revision/figures/"
ggsave(
  paste0(fig__path, "NoteS1FigA.png"),
  plot = figA, device = "png",
  width = 6, height = 6, dpi = 300)

