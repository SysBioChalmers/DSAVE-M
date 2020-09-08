#install.packages("jackstraw")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("lfa")
#BiocManager::install("qvalue")



library("jackstraw")
browseVignettes("jackstraw")



#the vignette

library(corpcor)
set.seed(1)
B = c(runif(100, min = 0.1, max = 1), rep(0, 900))
L = c(rep(1, 10), rep(-1, 10))
L = L/sd(L)
E = matrix(rnorm(1000 * 20), nrow = 1000)
Y = B %*% t(L) + E
dim(Y)

PA = permutationPA(Y, B = 10, threshold = 0.05)

plot(PA$p, pch = 20, main = "Permutation Parallel Analysis P-values",
     ylab = "P-values", xlab = "Principal Component")


svd.out = fast.svd(Y)
par(mfrow = c(2, 1))
plot(svd.out$d^2/sum(svd.out$d^2), pch = 20, main = "The scree plot",
     xlab = "PC", ylab = "Percent Variance Explained")
plot(svd.out$d[1] * svd.out$v[, 1], pch = 20, main = "1st PC",
     xlab = "Observation", ylab = "Magnitude")


