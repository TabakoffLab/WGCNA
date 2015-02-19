###################################################################
# Get Module Preservation Statistics Across All Other Brain Areas #
###################################################################

setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/acrossTissueCompare")
rm(list=ls())
options(stringsAsFactors = FALSE);
library(WGCNA)

###################
# Brain vs. Heart #
###################
load("heartExpr_RIstrainMeans.Rdata")
load("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/heritExprs_strainMeans.Rdata")

key = read.csv(file="/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/moduleMembership.csv", header=TRUE)
head(key)

labels_4analysis = key[,"moduleLabels"]
length(labels_4analysis)
table(labels_4analysis)
objects()

table(rownames(heritMeans)==rownames(heartMeans_RI))
table(rownames(heritMeans)==key$TC)

setLabels = c("Brain", "Heart");
#define your expression data (rows=strains, columns=transcriptClusters)
multiExpr = list(Brain = list(data=t(heritMeans)), Heart = list(data=t(heartMeans_RI)));
#define your module identification (based on brain);
multiColor = list(Brain = labels_4analysis);

system.time( {
mp = modulePreservation(multiExpr, multiColor,
networkType = "unsigned", corFnc = "bicor", referenceNetworks = 1,
nPermutations = 200,
randomSeed = 1,
maxModuleSize = 300,
maxGoldModuleSize = 300,
quickCor = 0,
verbose = 3)
} );

summary(mp)
summary(mp$preservation)
summary(mp$preservation$Z)
summary(mp$preservation$log.p$ref.Brain)
summary(mp$preservation$log.p$ref.Brain$inColumnsAlsoPresentIn.Heart)

Zsummary = mp$preservation$Z$ref.Brain$inColumnsAlsoPresentIn.Heart[,"Zsummary.pres"]
log.p = mp$preservation$log.p$ref.Brain$inColumnsAlsoPresentIn.Heart[,"log.psummary.pres"]
p.value = exp(log.p)

Zsummary.BrainRef.vHeart = cbind(Zsummary, log.p, p.value)
rownames(Zsummary.BrainRef.vHeart) = rownames(mp$preservation$Z$ref.Brain$inColumnsAlsoPresentIn.Heart)

save(mp, Zsummary.BrainRef.vHeart, file="Zsummary.BrainRef.vHeart.Rdata")
write.csv(Zsummary.BrainRef.vHeart, file="Zsummary.BrainRef.vHeart.csv")

