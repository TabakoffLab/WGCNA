#####################################################################
# Make Datasets Based on Brain RNAseq derived TC (20141218 version) #
# Do This for Heart, Liver and BAT 									#
#####################################################################


#Set up workspace;
rm(list=ls())
setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218")
library(WGCNA)
options(stringsAsFactors=FALSE)

#load in Brain (version 20141218) RNAseq derived TC information;
key = read.csv(file="PStoCollapsedTransciptClusterKey.csv")
length(unique(key[,"transcriptClusterCollapsed"]))
#but this is the key before the heritability filter....
herits = read.table(file="herits.clustersTry20141218.HXB.brain.txt",sep="\t",header=TRUE)

#only get key for heritable TC;
key_herits = merge(key, herits, by.x="transcriptClusterCollapsed", by.y="transcript")
dim(key_herits)
key_want = key_herits[which(key_herits[,"heritability"]>0.25),]
length(unique(key_want[,"transcriptClusterCollapsed"]))
#12,450 unique TC (62,103 probesets)

setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/acrossTissueCompare")


#########
# Heart #
#########
heartExpr = read.table(file="/data2/saba/ForPhenoGen/HXB.BXH.Heart.fullPS/rma.fullPS.HXB_BXH.heart.PhenoGen.txt",sep="\t",header=TRUE, row.names=1)

#get dataset with only the PS to collapse
heartWant = merge(as.matrix(key_want[,"probeset"]), heartExpr, by.x=1, by.y=0)
rownames(heartWant) = heartWant[,1]
heartWant = heartWant[,-1]

#collapse down into brain_RNAseq_TC;
heartCollapsed_byBrain20141218 = collapseRows(heartWant, key_want[,"transcriptClusterCollapsed"], key_want[,"probeset"], method="ME", thresholdCombine=NA)
length(heartCollapsed_byBrain20141218)
dim(heartCollapsed_byBrain20141218$datETcollapsed)
save(heartCollapsed_byBrain20141218, file="heartCollapsed_byBrain20141218.Rdata") 

#get strain means;
load("heartCollapsed_byBrain20141218.Rdata") 
exprHeart = heartCollapsed_byBrain20141218$datETcollapsed
strainHeart = sapply(strsplit(colnames(exprHeart), split="_", fixed=TRUE), "[[", 1)

heartMeans = t(apply(exprHeart, 1, function(a) lm(a~as.factor(strainHeart) -1)$coefficients))
colnames(heartMeans) = sapply(strsplit(colnames(heartMeans), split=")", fixed=TRUE), "[[", 2)

toDelete = c("BN.LX","PD","SHR.H","SHR.lj","SHR.Lx","WKY.Lj")
sig = c()
for(i in toDelete){
	sig = c(sig, which(colnames(heartMeans)==i))
}

length(sig)
heartMeans_RI = heartMeans[,-sig]
save(heartMeans_RI, file="heartExpr_RIstrainMeans.Rdata")


#########
# Liver #
#########
liverExpr = read.table(file="/data2/saba/ForPhenoGen/HXB.BXH.Liver.fullPS/rma.fullPS.HXB_BXH.liver.PhenoGen.txt",sep="\t",header=TRUE, row.names=1)

#get dataset with only the PS to collapse
liverWant = merge(as.matrix(key_want[,"probeset"]), liverExpr, by.x=1, by.y=0)
rownames(liverWant) = liverWant[,1]
liverWant = liverWant[,-1]

#collapse down into brain_RNAseq_TC;
liverCollapsed_byBrain20141218 = collapseRows(liverWant, key_want[,"transcriptClusterCollapsed"], key_want[,"probeset"], method="ME", thresholdCombine=NA)
length(liverCollapsed_byBrain20141218)
dim(liverCollapsed_byBrain20141218$datETcollapsed)
save(liverCollapsed_byBrain20141218, file="liverCollapsed_byBrain20141218.Rdata") 

#get strain means;
load("liverCollapsed_byBrain20141218.Rdata") 
exprliver = liverCollapsed_byBrain20141218$datETcollapsed
strainliver = sapply(strsplit(colnames(exprliver), split="_", fixed=TRUE), "[[", 1)

liverMeans = t(apply(exprliver, 1, function(a) lm(a~as.factor(strainliver) -1)$coefficients))
colnames(liverMeans) = sapply(strsplit(colnames(liverMeans), split=")", fixed=TRUE), "[[", 2)

toDelete = c("BN.LX","PD","SHR.H","SHR.lj","SHR.Lx","WKY.Lj")
sig = c()
for(i in toDelete){
	sig = c(sig, which(colnames(liverMeans)==i))
}

length(sig)
liverMeans_RI = liverMeans[,-sig]
save(liverMeans_RI, file="liverExpr_RIstrainMeans.Rdata")


#################
# Brown Adipose #
#################
BATExpr = read.table(file="/data2/saba/ForPhenoGen/HXB.BXH.BrownAdipose.fullPS/rma.fullPS.HXB_BXH.BAT.PhenoGen.txt",sep="\t",header=TRUE, row.names=1)

#get dataset with only the PS to collapse
BATWant = merge(as.matrix(key_want[,"probeset"]), BATExpr, by.x=1, by.y=0)
rownames(BATWant) = BATWant[,1]
BATWant = BATWant[,-1]

#collapse down into brain_RNAseq_TC;
BATCollapsed_byBrain20141218 = collapseRows(BATWant, key_want[,"transcriptClusterCollapsed"], key_want[,"probeset"], method="ME", thresholdCombine=NA)
length(BATCollapsed_byBrain20141218)
dim(BATCollapsed_byBrain20141218$datETcollapsed)
save(BATCollapsed_byBrain20141218, file="BATCollapsed_byBrain20141218.Rdata") 


#get strain means;
load("BATCollapsed_byBrain20141218.Rdata") 
exprBAT = BATCollapsed_byBrain20141218$datETcollapsed
strainBAT = sapply(strsplit(colnames(exprBAT), split="_", fixed=TRUE), "[[", 1)

BATMeans = t(apply(exprBAT, 1, function(a) lm(a~as.factor(strainBAT) -1)$coefficients))
colnames(BATMeans) = sapply(strsplit(colnames(BATMeans), split=")", fixed=TRUE), "[[", 2)

toDelete = c("BN.Lx","Pd","SHR.H","SHR.lj","SHR.Lx","WKY.lj")
sig = c()
for(i in toDelete){
	sig = c(sig, which(colnames(BATMeans)==i))
}

length(sig)
BATMeans_RI = BATMeans[,-sig]
save(BATMeans_RI, file="BATExpr_RIstrainMeans.Rdata")



