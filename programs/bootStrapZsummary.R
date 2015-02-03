##########################################################################
# Make Bootstrap Datasets for Brain RNAseq derived TC (20141218 version) #
# Get an idea on the variance of the Zsummary scores 			         #
##########################################################################

setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/acrossTissueCompare")
rm(list=ls())
options(stringsAsFactors = FALSE)
library(WGCNA)

load("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/heritExprs_strainMeans.Rdata")


##########################
# Get bootstrap datasets #
##########################
toPull = list()
for(i in 1:100){
	toPull[[i]] = sample(colnames(heritMeans), ncol(heritMeans), replace=TRUE)
}

brainBootStraps = list()
for(i in 1:10){
	brainBootStraps[[i]] = t(heritMeans[,toPull[[i]]])
	rownames(brainBootStraps[[i]]) = c(1:nrow(brainBootStraps[[i]]))
}

##############################################################
# Test Out How Long 10 Bootstrap Samples Would Take to Get a # 
# Z summary Score for 1 medium sized module 				 #
##############################################################

#load WGCNA results;
#key = read.csv(file="/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/moduleMembership.csv", header=TRUE)
key = read.csv(file="moduleMembership_withTestLabels.csv", header=TRUE)
head(key)
sizes = as.vector(table(key$moduleLabels))
summary(sizes)

#NEWlabel is 1 = module 36 (yellowgreen, size= 20) and 2 = module 37 (skyblue3, size=20);
sizes = as.vector(table(key$NEWlabel))
summary(sizes)

#Make sure we have the same order of genes in the datasets as well as in the label key;
table(colnames(t(heritMeans)) == key$TC)

getMP = function(x){
	setLabels = c("RefArea", "OtherArea");
	multiExpr = list(RefArea = list(data=t(heritMeans)), OtherArea = list(data=x));
	multiColor = list(RefArea = key$NEWlabel);

	mp = modulePreservation(multiExpr, multiColor,
		networkType = "unsigned", corFnc = "bicor", referenceNetworks = 1,
		nPermutations = 200,
		randomSeed = 1,
		maxModuleSize = 300,
		maxGoldModuleSize = 300,
		quickCor = 0,
		verbose = 3)

	return(mp)
}

MPinNiceFormat = function(x){
	Zsummary = x$preservation$Z$ref.RefArea$inColumnsAlsoPresentIn.OtherArea[,"Zsummary.pres"]
	log.p = x$preservation$log.p$ref.RefArea$inColumnsAlsoPresentIn.OtherArea[,"log.psummary.pres"]
	p.value = exp(log.p)

	Zsummary.Brain.vBoot = cbind(Zsummary, log.p, p.value)
	rownames(Zsummary.Brain.vBoot) = rownames(x$preservation$Z$ref.RefArea$inColumnsAlsoPresentIn.OtherArea)

	return(Zsummary.Brain.vBoot)	
}

system.time( {
boots1_to10_mp = lapply(brainBootStraps, getMP)
} );

#Took 2784.568 seconds to run;
2784.568/60
#46.4 minutes;

boots1_to10_mp_pretty = lapply(boots1_to10_mp, MPinNiceFormat)


