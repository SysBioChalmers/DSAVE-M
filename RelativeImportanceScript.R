#run install.packages("relaimpo") once before running this script
library(relaimpo)

#Set current directory to that of the source file; the path to the data is relative from there
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Read data
regrData = read.table(file="../TempData/RegrData.txt", header=T, sep="\t", row.names = 1)
head(regrData)
dim(regrData)

# Create linear regression model
fit_regrData = lm(score ~ ds_lc + ds_bc + ds_68k + ct_mon + ct_b + ct_nk + ct_mal + ti_tumlung + ti_tumbreast + ti_healtbreast + ti_lymphnode, data = regrData)
summary(fit_regrData)
# Run relative importance analysis
results = calc.relimp(fit_regrData, type = "lmg", rank = TRUE, diff = TRUE, rela = TRUE)$lmg
#Look at results convert to %. Copy these values to excel.
results*100
# Plot (for an overview only, the graphs are created in Excel)
barplot(sort(results*100, decreasing = TRUE), las=2, ylim=c(0,100), ylab="Relative explained variance (%)")


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
results2*100
# Plot (for an overview only, the graphs are created in Excel)
barplot(sort(results2*100, decreasing = TRUE), las=2, ylim=c(0,100), ylab="Relative exaplained variance (%)")
