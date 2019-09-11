#run install.packages("relaimpo") once before running this script
library(relaimpo)

#Set current directory to that of the source file; the path to the data is relative from there
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Read data
regrData = read.table(file="../../TempData/RegrData.txt", header=T, sep="\t", row.names = 1)
head(regrData)
dim(regrData)

# Create linear regression model
fit_regrData = lm(score ~ ds_lc + ds_bc + ds_68k + ct_mon + ct_b + ct_nk + ct_mal + ti_tumlung + ti_tumbreast + ti_healtbreast + ti_lymphnode, data = regrData)
summary(fit_regrData)
# Run relative importance analysis
results = calc.relimp(fit_regrData, type = "lmg", rank = TRUE, diff = TRUE, rela = TRUE)$lmg
#Look at results convert to %. Copy these values to excel.
#results*100
# Plot (for an overview only, the graphs are created in Excel)
#barplot(sort(results*100, decreasing = TRUE), las=2, ylim=c(0,100), ylab="Relative explained variance (%)")

#scale the results to balance out the unequal number of samples per factor
#These are the values used in the paper!
sampPerFactor = colSums(regrData)
#remove score and intercept
sampPerFactor = sampPerFactor[-1]
sampPerFactor = sampPerFactor[-1]
#divide results with samples per factor and renormalize
scaledResults = results / sampPerFactor
scaledResults = scaledResults / sum(scaledResults)
scaledResults*100
barplot(sort(scaledResults*100, decreasing = TRUE), las=2, ylim=c(0,100), ylab="Relative explained variance (%)")

################### Remove breast cancer and run again
rows = row.names(regrData)
sel = substr(rows,1,2) == "bc"
regrData2 = regrData[!sel, ]
head(regrData2)
dim(regrData2)

# Create linear regression model
fit_regrData2 = lm(score ~ ds_lc + ds_68k + ct_mon + ct_b + ct_nk + ct_mal + ti_tumlung, data = regrData2)
summary(fit_regrData2)
# Run relative importance analysis
results2 = calc.relimp(fit_regrData2, type = "lmg", rank = TRUE, diff = TRUE, rela = TRUE)$lmg
#Look at results
#results2*100
# Plot (for an overview only, the graphs are created in Excel)
#barplot(sort(results2*100, decreasing = TRUE), las=2, ylim=c(0,100), ylab="Relative exaplained variance (%)")

#Scale the results to balance out the unequal number of samples per factor
#These are the values used in the paper!
sampPerFactor2 = sampPerFactor[c(1,3,4,5,6,7,8)]
#divide results with samples per factor and renormalize
scaledResults2 = results2 / sampPerFactor2
scaledResults2 = scaledResults2 / sum(scaledResults2)
scaledResults2*100
barplot(sort(scaledResults2*100, decreasing = TRUE), las=2, ylim=c(0,100), ylab="Relative explained variance (%)")

