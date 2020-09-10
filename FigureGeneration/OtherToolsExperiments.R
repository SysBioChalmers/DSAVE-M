#THis file contains a lot of experimental code which was not in the end used for any figure generation

#install.packages("jackstraw")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("lfa")
#BiocManager::install("qvalue")



library(jackstraw)
library(ggplot2)

#dat <- read.csv("ncomms14049-s2.csv")
impFolder = "C:/Work/MatlabCode/components/DSAVE-M/ImportableData/"
dat <- read.table(paste0(impFolder,"41467_2017_BFncomms14049_MOESM829_ESM.csv"), row.names = NULL, sep=";", header = T)

dat <- dat[dat$Sample == "50%:50% Jurkat:293T", ]
dat <- dat[, -1]
Jurkat293T <- dat

set.seed(0)
dat.plot <- dat
k <- kmeans(dat, 2, iter.max = 100, nstart = 100)
dat.plot$k <- k$cluster
g <- ggplot(dat.plot, aes(PC1, PC2, col = as.factor(k))) + geom_point(size = 0.2,
                                                                      alpha = 1) + theme_classic() + scale_color_manual(values = c("navyblue", "darkorange"))
print(g)



set.seed(0)
js <- jackstraw_kmeans(as.matrix(dat), kmeans.dat = k, s = round(nrow(dat) * 0.1),
                       B = 1000, verbose = FALSE, pool = FALSE)

dat.plot <- dat
dat.plot$k <- k$cluster
dat.plot$Pvalue <- js$p.F

ggplot(dat.plot, aes(Pvalue)) + geom_histogram(bins = 30) + theme_minimal()


g.pvalue <- ggplot(dat.plot, aes(PC1, Pvalue)) + geom_point(aes(col = as.factor(k),
                                                                shape = as.factor(k)), size = 0.6) + scale_colour_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                              order = 1), labels = c("1", "2"), values = c("navyblue", "darkorange")) + scale_shape_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                                                                                                                order = 1), labels = c("1", "2"), values = c(16, 17)) + theme_minimal() + theme(legend.position = "none")
print(g.pvalue)

library(qvalue)
# lambda estimated from fitting a cubic smoothing spline
qvalue::pi0est(js$p.F)$pi0

library(mutoss)
mutoss::BR_pi0_est(js$p.F, alpha = 0.3)$pi0

pi0_avg = (qvalue::pi0est(js$p.F)$pi0 + mutoss::BR_pi0_est(js$p.F, alpha = 0.3)$pi0)/2

pipout = pip(js$p.F, group = k$cluster, pi0 = pi0_avg)

library(Matrix)
require(gridExtra)
library(cowplot)

dat.plot <- dat
dat.plot$k <- k$cluster
dat.plot$PIP <- pipout
dat.plot$Pvalue <- js$p.F

g.pip <- ggplot(dat.plot, aes(PC1, PIP)) + geom_point(aes(col = as.factor(k), shape = as.factor(k)),
                                                      size = 0.6) + scale_colour_manual(guide = guide_legend(title = "Cluster", order = 1),
                                                                                        labels = c("1", "2"), values = c("navyblue", "darkorange")) + scale_shape_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                                              order = 1), labels = c("1", "2"), values = c(16, 17)) + theme_minimal() + theme(legend.position = "none")
g.pvalue.pip <- ggplot(dat.plot, aes(Pvalue, PIP)) + geom_point(aes(col = as.factor(k),
                                                                    shape = as.factor(k)), size = 0.6) + geom_hline(yintercept = 0.8, col = "red",
                                                                                                                    lty = "dashed") + scale_x_continuous(limits = c(0, 1)) + scale_colour_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                                                                      order = 1), labels = c("1", "2"), values = c("navyblue", "darkorange")) + scale_shape_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                                                                                                                                                                                        order = 1), labels = c("1", "2"), values = c(16, 17)) + theme_minimal() + theme(legend.position = "none")
print(g.pvalue.pip)

require(scales)
highcol <- "black"
lowcol <- "white"

g.pca <- ggplot(dat.plot, aes(PC1, PC2, alpha = PIP)) + geom_point(aes(alpha = PIP,
                                                                       col = as.factor(k), shape = as.factor(k)), size = 0.6) + scale_colour_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                         order = 1, override.aes = list(size = 1.5)), labels = c("1", "2"), values = c("navyblue",
                                                                                                                                                                                                                                                       "darkorange")) + scale_shape_manual(guide = guide_legend(title = "Cluster", order = 1),
                                                                                                                                                                                                                                                                                           labels = c("1", "2"), values = c(16, 17)) + geom_point(aes(fill = PIP), alpha = 0) +
  scale_fill_gradient(high = highcol, low = lowcol) + guides(alpha = F) + theme_minimal() +
  theme(legend.position = "right")
print(g.pca)


#test to show the ones with PIP < 0.8
dat.plot2 = dat.plot
dat.plot2$filt = dat.plot$PIP < 0.8

g.pca2 <- ggplot(dat.plot2, aes(PC1, PC2, alpha = PIP)) + geom_point(aes(alpha = PIP,
                                                                       col = filt, shape = as.factor(k)), size = 0.6) + scale_colour_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                         order = 1, override.aes = list(size = 1.5)), labels = c("1", "2"), values = c("navyblue",
                                                                                                                                                                                                                                                       "darkorange")) + scale_shape_manual(guide = guide_legend(title = "Cluster", order = 1),
                                                                                                                                                                                                                                                                                           labels = c("1", "2"), values = c(16, 17)) + geom_point(aes(fill = PIP), alpha = 0) +
  scale_fill_gradient(high = highcol, low = lowcol) + guides(alpha = F) + theme_minimal() +
  theme(legend.position = "right")
print(g.pca2)


#now, the seurat example
#library(ggplot2)
#library(DSAVE)
library(progress)
library(Seurat)
library("plotly")
library("downloader")

set.seed(0)
extrDir <- downloadData("https://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz", "PBMC68k")
dataDir = paste0(extrDir,"/filtered_matrices_mex/hg19")
pbmc <- Read10X(data.dir = dataDir)
#only use the 20000 first cells to speed up the demo

pbmc = pbmc[,1:20000]
#pbmcCellTypes = ctPbmc68k[1:20000]

#further filter out a fraction of the cells from cell types we will not use to speed up the demo
#toRem = (pbmcCellTypes == "CD8+/CD45RA+ Naive Cytotoxic") | (pbmcCellTypes == "CD8+ Cytotoxic T")
#toRem[15000:20000] = F
#pbmc = pbmc[,!toRem]
#pbmcCellTypes = pbmcCellTypes[!toRem]

#run through Seurat
so = CreateSeuratObject(counts = pbmc, project = "pbmc", min.cells = 0, min.features = 0)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)

#Finds the most variable genes
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)

#Make mean and variance the same for all genes:
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)

#Principal component analysis
so <- RunPCA(so, features = VariableFeatures(object = so))

#so[["extCellTypes"]] = pbmcCellTypes

so <- RunUMAP(so, dims = 1:10)
so <- FindNeighbors(so, dims = 1:10)
so <- FindClusters(so, resolution = 0.5)
DimPlot(so, reduction = "umap")
#DimPlot(so, reduction = "umap", group.by = "extCellTypes")
c7 = extractSeuratData(so, 7)
#also investigate the clustering from Seurat directly
set.seed(0)
c7_divergence <- DSAVEGetSingleCellDivergence(c7,
                                                  minUMIsPerCell = 200, tpmLowerBound = 0,
                                                  iterations = 15,silent=FALSE)

DSAVEPlotDivergence(c7, c7_divergence)

#now test to use kmeans clustering since this is easier to use with Jackstraw:
set.seed(0)

#first filter empty genes to reduce the data
#totCounts = sparse_Sums(pbmc, rowSums = T)
#sum(totCounts == 0)

#extract the PCs from Seurat
so@reductions
pcaData = Embeddings(so, reduction = "pca")[, 1:10]

k <- kmeans(pcaData, 10, iter.max = 100, nstart = 100)
dat.plot = as.data.frame(pcaData)
#dat.plot$k <- as.factor(as.numeric(k$cluster))
g <- ggplot(dat.plot, aes(PC_1, PC_2, col = k)) + geom_point(size = 0.2,
                                                                      alpha = 1) + theme_classic()
print(g)

so[["kmeansClusters"]] = k$cluster
DimPlot(so, reduction = "umap", group.by = "kmeansClusters")


pbmc_js <- jackstraw_kmeans(pcaData, kmeans.dat = k, s = round(nrow(pcaData) * 0.1),
                       B = 1000, verbose = FALSE, pool = FALSE)

dat.plot$k <- k$cluster
dat.plot$Pvalue <- pbmc_js$p.F

ggplot(dat.plot, aes(Pvalue)) + geom_histogram(bins = 30) + theme_minimal()

g.pvalue <- ggplot(dat.plot, aes(PC_1, Pvalue)) + geom_point(aes(col = as.factor(k)), size = 0.6)
print(g.pvalue)

soSub = so[, so[["kmeansClusters"]] == 8]
DimPlot(soSub, reduction = "umap", group.by = "kmeansClusters")

## Investigate divergence on cluster 8
#########################
datsub = pbmc[, as.vector(unname(so[["kmeansClusters"]] == 8))]

set.seed(0)
sub_divergence <- DSAVEGetSingleCellDivergence(datsub,
                                              minUMIsPerCell = 200, tpmLowerBound = 0,
                                              iterations = 15)

DSAVEPlotDivergence(datsub, sub_divergence)
#191,292
unwhich <- function (which, dim = max(which)) { #copied from https://stat.ethz.ch/pipermail/r-help/2001-July/013815.html
  y <- array(logical(length(which)), dim = dim)
  y[which] <- TRUE
  y
}

#GNLY, GZMB, NKG7 (cytotoxic T cells or NK cells)
soSub[["markedCells"]] = unwhich(c(74,191,292,418,422,446), dim(datsub)[2])

soSub["GZMB", 292]$nCount_RNA # test, looks ok

DimPlot(soSub, reduction = "umap", group.by = "markedCells")
#so, our two cells with NK cell markers are in the middle of the cluster (which is not NK cells)

#Now test with Jackstraw


library(qvalue)
# lambda estimated from fitting a cubic smoothing spline
qvalue::pi0est(pbmc_js$p.F)$pi0

library(mutoss)
mutoss::BR_pi0_est(pbmc_js$p.F, alpha = 0.3)$pi0

pi0_avg = (qvalue::pi0est(pbmc_js$p.F)$pi0 + mutoss::BR_pi0_est(pbmc_js$p.F, alpha = 0.3)$pi0)/2

pipout = pip(pbmc_js$p.F, group = k$cluster, pi0 = pi0_avg)

pipSub = pipout[as.vector(unname(so[["kmeansClusters"]] == 8))]
pipSub[418]

soSub[["lowPip"]] = pipSub < 0.8
sum(pipSub < 0.8)
which(pipSub < 0.8)

DimPlot(soSub, reduction = "umap", group.by = "lowPip")


library(Matrix)
require(gridExtra)
library(cowplot)

dat.plot <- dat
dat.plot$k <- k$cluster
dat.plot$PIP <- pipout
dat.plot$Pvalue <- js$p.F

g.pip <- ggplot(dat.plot, aes(PC1, PIP)) + geom_point(aes(col = as.factor(k), shape = as.factor(k)),
                                                      size = 0.6) + scale_colour_manual(guide = guide_legend(title = "Cluster", order = 1),
                                                                                        labels = c("1", "2"), values = c("navyblue", "darkorange")) + scale_shape_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                                              order = 1), labels = c("1", "2"), values = c(16, 17)) + theme_minimal() + theme(legend.position = "none")
g.pvalue.pip <- ggplot(dat.plot, aes(Pvalue, PIP)) + geom_point(aes(col = as.factor(k),
                                                                    shape = as.factor(k)), size = 0.6) + geom_hline(yintercept = 0.8, col = "red",
                                                                                                                    lty = "dashed") + scale_x_continuous(limits = c(0, 1)) + scale_colour_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                                                                      order = 1), labels = c("1", "2"), values = c("navyblue", "darkorange")) + scale_shape_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                                                                                                                                                                                        order = 1), labels = c("1", "2"), values = c(16, 17)) + theme_minimal() + theme(legend.position = "none")
print(g.pvalue.pip)

require(scales)
highcol <- "black"
lowcol <- "white"

g.pca <- ggplot(dat.plot, aes(PC1, PC2, alpha = PIP)) + geom_point(aes(alpha = PIP,
                                                                       col = as.factor(k), shape = as.factor(k)), size = 0.6) + scale_colour_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                         order = 1, override.aes = list(size = 1.5)), labels = c("1", "2"), values = c("navyblue",
                                                                                                                                                                                                                                                       "darkorange")) + scale_shape_manual(guide = guide_legend(title = "Cluster", order = 1),
                                                                                                                                                                                                                                                                                           labels = c("1", "2"), values = c(16, 17)) + geom_point(aes(fill = PIP), alpha = 0) +
  scale_fill_gradient(high = highcol, low = lowcol) + guides(alpha = F) + theme_minimal() +
  theme(legend.position = "right")
print(g.pca)


#test to show the ones with PIP < 0.8
dat.plot2 = dat.plot
dat.plot2$filt = dat.plot$PIP < 0.8

g.pca2 <- ggplot(dat.plot2, aes(PC1, PC2, alpha = PIP)) + geom_point(aes(alpha = PIP,
                                                                         col = filt, shape = as.factor(k)), size = 0.6) + scale_colour_manual(guide = guide_legend(title = "Cluster",
                                                                                                                                                                   order = 1, override.aes = list(size = 1.5)), labels = c("1", "2"), values = c("navyblue",
                                                                                                                                                                                                                                                 "darkorange")) + scale_shape_manual(guide = guide_legend(title = "Cluster", order = 1),
                                                                                                                                                                                                                                                                                     labels = c("1", "2"), values = c(16, 17)) + geom_point(aes(fill = PIP), alpha = 0) +
  scale_fill_gradient(high = highcol, low = lowcol) + guides(alpha = F) + theme_minimal() +
  theme(legend.position = "right")
print(g.pca2)


#now test scReclassify
#devtools::install_github("SydneyBioX/scReClassify", build_opts = c("--no-resave-data", "--no-manual"))
#install.packages("mclust")
library(scReClassify)

pcaData2 = t(Embeddings(so, reduction = "pca")[,1:20])

library(tictoc)
tic()
reclassifyData <- multiAdaSampling(pcaData2, k$cluster, seed = 1, classifier = "svm", percent = 1, L = 10)
toc()

#extract the cluster 8 cells
reclassified = reclassifyData$final[k$cluster == 8]
idsReclassified = which(reclassified != 8)


# check examples
idx <- which(cellTypes.reclassify$final != cellTypes)
library(dplyr)
cbind(original=cellTypes[idx], reclassify=cellTypes.reclassify$final[idx]) %>%
  DT::datatable()

c1 <- dat.processed[, which(cellTypes=="Endothelial Cell")]
c2 <- dat.processed[, which(cellTypes=="Erythrocyte")]
c3 <- dat.processed[, which(cellTypes=="Hepatoblast")]
c4 <- dat.processed[, which(cellTypes=="Macrophage")]
c5 <- dat.processed[, which(cellTypes=="Megakaryocyte")]
c6 <- dat.processed[, which(cellTypes=="Mesenchymal Cell")]
cs <- rainbow(length(table(cellTypes)))

# (example 1 E13.5_C20)
#####
par(mfrow=c(1,2))
marker <- End[1]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker, names=c("Reclassified", "Orignal", "Others", "Others", "Others", "Others"), las=2)
points(1, dat.processed[marker, which(colnames(dat.processed) %in% "E13.5_C20")], pch=16, col="red", cex=2)

marker <- End[2]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker, names=c("Reclassified", "Orignal", "Others", "Others", "Others", "Others"), las=2)
points(1, dat.processed[marker, which(colnames(dat.processed) %in% "E13.5_C20")], pch=16, col="red", cex=2)















# Verification by marker genes
End <- c("KDR", "LYVE1")

# check examples
idx <- which(cellTypes.reclassify$final != cellTypes)






data(GSE87795_liver.development.data)
dat <- GSE87795_liver.development.data$data
cellTypes <- GSE87795_liver.development.data$cellTypes

# number of clusters
nCs <- length(table(cellTypes))

# This demo dataset is already pre-processed
dat.processed = dat
dat.selected = matPCs(dat.processed, 0.7)

lab <- cellTypes

set.seed(1)
noisyCls <- function(dat, rho, cls.truth){
  cls.noisy <- cls.truth
  names(cls.noisy) <- colnames(dat)
  for(i in 1:length(table(cls.noisy))) {
    # class label starts from 0
    if (i != length(table(cls.noisy))) {
      cls.noisy[sample(which(cls.truth == names(table(cls.noisy))[i]), floor(sum(cls.truth == names(table(cls.noisy))[i]) * rho))] <- names(table(cls.noisy))[i+1]
    } else {
      cls.noisy[sample(which(cls.truth == names(table(cls.noisy))[i]), floor(sum(cls.truth == names(table(cls.noisy))[i]) * rho))] <- names(table(cls.noisy))[1]
    }
  }

  print(sum(cls.truth != cls.noisy))
  return(cls.noisy)
}

cls.noisy01 <- noisyCls(dat.selected, rho=0.1, lab)
cls.noisy02 <- noisyCls(dat.selected, rho=0.2, lab)
cls.noisy03 <- noisyCls(dat.selected, rho=0.3, lab)
cls.noisy04 <- noisyCls(dat.selected, rho=0.4, lab)
cls.noisy05 <- noisyCls(dat.selected, rho=0.5, lab)


###################################
# SVM
###################################
acc01 <- acc02 <- acc03 <- acc04 <- acc05 <- c()
ari01 <- ari02 <- ari03 <- ari04 <- ari05 <- c()
base <- "svm"

for(j in 1:10) {
  final <- multiAdaSampling(dat.selected, cls.noisy01, seed=j, classifier=base, percent=1, L=10)$final
  ari01 <- c(ari01, mclust::adjustedRandIndex(lab, final))
  acc01 <- c(acc01, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy02, seed=j, classifier=base, percent=1, L=10)$final
  ari02 <- c(ari02, mclust::adjustedRandIndex(lab, final))
  acc02 <- c(acc02, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy03, seed=j, classifier=base, percent=1, L=10)$final
  ari03 <- c(ari03, mclust::adjustedRandIndex(lab, final))
  acc03 <- c(acc03, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy04, seed=j, classifier=base, percent=1, L=10)$final
  ari04 <- c(ari04, mclust::adjustedRandIndex(lab, final))
  acc04 <- c(acc04, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy05, seed=j, classifier=base, percent=1, L=10)$final
  ari05 <- c(ari05, mclust::adjustedRandIndex(lab, final))
  acc05 <- c(acc05, bAccuracy(lab, final))
}

result = list(
  acc01 = acc01,
  acc02 = acc02,
  acc03 = acc03,
  acc04 = acc04,
  acc05 = acc05,
  ari01 = ari01,
  ari02 = ari02,
  ari03 = ari03,
  ari04 = ari04,
  ari05 = ari05
)

plot.new()
par(mfrow = c(1,2))
boxplot(acc01, acc02, acc03, acc04, acc05, col="lightblue", main="SVM Acc", ylim=c(0.45, 1))
points(x=1:5, y=c(bAccuracy(lab, cls.noisy01), bAccuracy(lab, cls.noisy02),
                  bAccuracy(lab, cls.noisy03), bAccuracy(lab, cls.noisy04),
                  bAccuracy(lab, cls.noisy05)), col="red3", pch=c(2,3,4,5,6), cex=1)
boxplot(ari01, ari02, ari03, ari04, ari05, col="lightblue", main="SVM ARI", ylim=c(0.25, 1))
points(x=1:5, y=c(mclust::adjustedRandIndex(lab, cls.noisy01), mclust::adjustedRandIndex(lab, cls.noisy02),
                  mclust::adjustedRandIndex(lab, cls.noisy03), mclust::adjustedRandIndex(lab, cls.noisy04),
                  mclust::adjustedRandIndex(lab, cls.noisy05)), col="red3", pch=c(2,3,4,5,6), cex=1)


# PCA procedure
dat.pc <- matPCs(dat.processed, 0.7)
dim(dat.pc)

# run scReClassify
cellTypes.reclassify <- multiAdaSampling(dat.pc, cellTypes, seed = 1, classifier = "svm", percent = 1, L = 10)

# Verification by marker genes
End <- c("KDR", "LYVE1")

# check examples
idx <- which(cellTypes.reclassify$final != cellTypes)
library(dplyr)
cbind(original=cellTypes[idx], reclassify=cellTypes.reclassify$final[idx]) %>%
  DT::datatable()

c1 <- dat.processed[, which(cellTypes=="Endothelial Cell")]
c2 <- dat.processed[, which(cellTypes=="Erythrocyte")]
c3 <- dat.processed[, which(cellTypes=="Hepatoblast")]
c4 <- dat.processed[, which(cellTypes=="Macrophage")]
c5 <- dat.processed[, which(cellTypes=="Megakaryocyte")]
c6 <- dat.processed[, which(cellTypes=="Mesenchymal Cell")]
cs <- rainbow(length(table(cellTypes)))

# (example 1 E13.5_C20)
#####
par(mfrow=c(1,2))
marker <- End[1]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker, names=c("Reclassified", "Orignal", "Others", "Others", "Others", "Others"), las=2)
points(1, dat.processed[marker, which(colnames(dat.processed) %in% "E13.5_C20")], pch=16, col="red", cex=2)

marker <- End[2]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker, names=c("Reclassified", "Orignal", "Others", "Others", "Others", "Others"), las=2)
points(1, dat.processed[marker, which(colnames(dat.processed) %in% "E13.5_C20")], pch=16, col="red", cex=2)
