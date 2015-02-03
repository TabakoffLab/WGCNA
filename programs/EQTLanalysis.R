rm(list=ls())
setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218")
library(qtl)

##Load in data
MEs <- read.cross("csv", "", "forAnalysis_MEs.csv", genotypes = c("H", "B"))
summary(MEs)


#########################
# QTL marker regression #
#########################
nphe=(nphe(MEs)-2)
out_MEs.mr <- scanone(MEs, pheno.col=c(1:nphe), method="mr")
max(out_MEs.mr)

dim(out_MEs.mr)
out_MEs.mr[1:5, 1:5]
head(out_MEs.mr)


maxLODs <- apply(out_MEs.mr[,3:ncol(out_MEs.mr)],2,function(a) out_MEs.mr[which.max(a),1:2])
results <- data.frame(probeset_id = names(maxLODs),snp_id = unlist(lapply(maxLODs,function(a) rownames(a))),
chr = unlist(lapply(maxLODs,function(a) a[,1])),pos = unlist(lapply(maxLODs,function(a) a[,2])),maxLOD = apply(out_MEs.mr[,3:ncol(out_MEs.mr)],2,max))

################
# get p-values #
################
set.seed(12)
operms_corrMEs <- scanone(MEs, pheno.col=c(1:nphe), method="mr", n.perm=1000)	
dim(operms_corrMEs)
summary(operms_corrMEs)

pvals <- c()
for(i in 1:nrow(results)){
	pvalue <- sum(results[i, "maxLOD"] < operms_corrMEs[,i])/nrow(operms_corrMEs)
	pvals <- c(pvals, pvalue)
}

results$pvals <- pvals

save(out_MEs.mr, operms_corrMEs, results, file="EQTLs.Rdata")
write.csv(results, file="HXB.RNAseqTC.EQTLs.20141218results.csv", row.names=FALSE) 
