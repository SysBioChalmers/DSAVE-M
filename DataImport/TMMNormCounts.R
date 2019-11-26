
library("edgeR")

folder = "C:/Work/MatlabCode/components/DSAVE-M/ImportableData" # to be replaced by each user
bulk_counts = read.table(file=paste0(folder, "/countsMatrix.txt"), header=T, sep="\t")

#sel_bulk_counts = bulk_counts[,32:39]
#sel_bulk_counts = bulk_counts[,c(30:32,34:37,39)]
sel_bulk_counts = bulk_counts[,c(30:36,39)]

#normalize using TMM
normFactors <-  calcNormFactors(sel_bulk_counts)

libSizes = colSums(sel_bulk_counts)
effectiveLibSizes = libSizes * normFactors


#I need to transpose the data matrix back and forth to get the
#row wise division to work...
scaledTMMs = t(t(sel_bulk_counts) * (10^6/ effectiveLibSizes))


#There are some issues with roundoff or similar making the mean of the sum of all 
#genes not to be exactly 10^6. Fix this:
sumPerSamp = colSums(scaledTMMs)
meanSampSum = mean(sumPerSamp)

scaledTMMs = scaledTMMs * 10^6/meanSampSum


#Test
colSums(scaledTMMs)
mean(colSums(scaledTMMs))

#Export
write.table(scaledTMMs, file=paste0(folder, "/scaledTMMMatrix.txt"), row.names=T, sep="\t")
