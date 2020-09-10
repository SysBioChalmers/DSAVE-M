#Evaluates DSAVE against two other tools: Jackstraw and scReClassify
#The DSAVE package is expected to be loaded

#install.packages("jackstraw")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("lfa")
#BiocManager::install("qvalue")



library(jackstraw)
library(ggplot2)


library(qvalue)
library(mutoss)
library(Matrix)
require(gridExtra)
library(cowplot)
require(scales)


#now, the seurat example

#library(DSAVE)
library(progress)
library(Seurat)
library("plotly")
library("downloader")
library(ggpubr)


fig__path = "Z:/projects/Partitioning of Cell-to-Cell Variation in Single-Cell RNA-Seq Data/paper/PLOS ONE/Revision/figures/"


set.seed(0)
extrDir <- downloadData("https://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz", "PBMC68k")
dataDir = paste0(extrDir,"/filtered_matrices_mex/hg19")
pbmc <- Read10X(data.dir = dataDir)

#only use the 20000 first cells to speed up the analysis
pbmc = pbmc[,1:20000]

#run through Seurat
so = CreateSeuratObject(counts = pbmc, project = "pbmc", min.cells = 0, min.features = 0)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)
so <- RunPCA(so, features = VariableFeatures(object = so))
so <- RunUMAP(so, dims = 1:10)
DimPlot(so, reduction = "umap")

#Use kmeans clustering instead of the built-in in Seurat clustering since this is easier to use with Jackstraw:
set.seed(0)

#first filter empty genes to reduce the data
#totCounts = sparse_Sums(pbmc, rowSums = T)
#sum(totCounts == 0)

#extract the PCs from Seurat
pcaData = Embeddings(so, reduction = "pca")[, 1:10]

k <- kmeans(pcaData, 10, iter.max = 100, nstart = 100)

so[["kmeansClusters"]] = k$cluster
figA1 = DimPlot(so, reduction = "umap", group.by = "kmeansClusters") + ggtitle("All clusters")
figA1

pbmc_js <- jackstraw_kmeans(pcaData, kmeans.dat = k, s = round(nrow(pcaData) * 0.1),
                       B = 1000, verbose = FALSE, pool = FALSE)

soSub = so[, so[["kmeansClusters"]] == 4]
figA2 = DimPlot(soSub, reduction = "umap") + theme(legend.position = "none") + ggtitle("Cluster 4")
figA2

figA = ggarrange(figA1, figA2, nrow=1, ncol=2,labels=c("I","II"))
figA

ggsave(
  paste0(fig__path, "NoteS3FigA.png"),
  plot = figA, device = "png",
  width = 7, height = 3, dpi = 300)



figB = FeaturePlot(so, c("NKG7", "LYZ", "CD14", "FCGR3A", "MS4A7", "FCER1A", "CST3"))
figB
ggsave(
  paste0(fig__path, "NoteS3FigB.png"),
  plot = figB, device = "png",
  width = 8, height = 8, dpi = 300)

#suggests FCGR3A+ Mono

## Investigate divergence on cluster 4
#########################
datsub = pbmc[, as.vector(unname(so[["kmeansClusters"]] == 4))]
dim(datsub)

set.seed(0)
sub_divergence <- DSAVEGetSingleCellDivergence(datsub,
                                              minUMIsPerCell = 200, tpmLowerBound = 0,
                                              iterations = 15)

#interactive plot for finding suitable cells
DSAVEPlotDivergence(datsub, sub_divergence)

unwhich <- function (which, dim = max(which)) { #copied from https://stat.ethz.ch/pipermail/r-help/2001-July/013815.html
  y <- array(logical(length(which)), dim = dim)
  y[which] <- TRUE
  y
}

#74,97,191,292,418,422, 446
#GNLY, GZMB, NKG7 (cytotoxic T cells or NK cells)
selectedCells = c(74,97,191,292,418,422,446)
soSub[["markedCells"]] = unwhich(selectedCells, dim(datsub)[2])

soSub["GZMB", 292]$nCount_RNA # test, looks ok

DimPlot(soSub, reduction = "umap", group.by = "markedCells")
#so, our two cells with NK cell markers are in the middle of the cluster (which is not NK cells)

#Now test with Jackstraw

# lambda estimated from fitting a cubic smoothing spline
qvalue::pi0est(pbmc_js$p.F)$pi0
mutoss::BR_pi0_est(pbmc_js$p.F, alpha = 0.3)$pi0
pi0_avg = (qvalue::pi0est(pbmc_js$p.F)$pi0 + mutoss::BR_pi0_est(pbmc_js$p.F, alpha = 0.3)$pi0)/2
pipout = pip(pbmc_js$p.F, group = k$cluster, pi0 = pi0_avg)
pipSub = pipout[as.vector(unname(so[["kmeansClusters"]] == 4))]
pipSub[418]

soSub[["lowPip"]] = pipSub < 0.8
sum(pipSub < 0.8)
divCellsJackstraw = which(pipSub < 0.8)
divCellsJackstraw
intersect(divCellsJackstraw, selectedCells)

DimPlot(soSub, reduction = "umap", group.by = "lowPip")


#now test scReclassify
#devtools::install_github("SydneyBioX/scReClassify", build_opts = c("--no-resave-data", "--no-manual"))
#install.packages("mclust")
library(scReClassify)

pcaData2 = t(Embeddings(so, reduction = "pca")[,1:20])

library(tictoc)
tic()
reclassifyData <- multiAdaSampling(pcaData2, k$cluster, seed = 1, classifier = "svm", percent = 1, L = 10)
toc()

#extract the cluster 4 cells
reclassified = reclassifyData$final[k$cluster == 4]
idsReclassified = which(reclassified != 4)
length(allIdsReclassified)
allIdsReclassified = which(reclassifyData$final != k$cluster)
hist(k$cluster[reclassifyData$final != k$cluster])
length(allIdsReclassified)


#Create figure C for S3
countsPerCell = textTinyR::sparse_Sums(datsub, rowSums=F)
dfC1 = data.frame(x=sub_divergence$divs, y=countsPerCell, NK = unwhich(selectedCells, dim(datsub)[2]))
sb = rep(20,dim(datsub)[2])
sb[dfC1$NK] = 6
#sb2 = as.factor(sb)
Divergent = factor(sb, labels = c("NK", "Other"))
dfC1b = dfC1[divCellsJackstraw,]

figC1 = ggplot(data=dfC1, aes(x=x, y=y, color=Divergent, shape=Divergent)) +
  geom_point() +
  geom_point(data=dfC1b, color="#000000", shape=1, size=3) +
  labs(title="Divergence of cluster 4", x="Divergence", y="UMI counts") +
  theme_bw()

figC1

UMAPCoords = t(Embeddings(soSub, reduction = "umap"))
dfC2 = data.frame(x=as.numeric(UMAPCoords[1,]), y=as.numeric(UMAPCoords[2,]))
dfC2b = dfC2[divCellsJackstraw,]

figC2 = ggplot(data=dfC2, aes(x=x, y=y, color=Divergent, shape=Divergent)) +
  geom_point() +
  geom_point(data=dfC2b, color="#000000", shape=1, size=3) +
  labs(title="Cluster 4, UMAP space", x="UMAP 1", y="UMAP 2") +
  theme_bw()
figC2

figC = ggarrange(figC1, figC2, nrow=1, ncol=2,labels=c("I","II"))
figC

ggsave(
  paste0(fig__path, "NoteS3FigC.png"),
  plot = figC, device = "png",
  width = 9, height = 4, dpi = 300)

