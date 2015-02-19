##Get a ME file with the Strains as rownames;
setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218")
rm(list=ls())
library(WGCNA)
options(strainsAsFactors=FALSE)

#load in MEs;
load("HXBbrain.RNAseqTC.WGCNAcontruct.20141218.Rdata")

#same strains as brain...
key = read.csv(file="/data2/kiemele/WGCNA/HXBheart_forPhenogen/20141231/origStrainToStrainKey.csv")

#get new ME file;
#attach the original strain names (becuase it is based on this order);
rownames(MEs) = key[,"origStrain"]

MEs_new = merge(key, MEs, by.x="origStrain", by.y=0)
rownames(MEs_new) = MEs_new[,"strain"]
MEs_new = MEs_new[order(MEs_new[,"strain"]),]

#just keep the rownames as an identifier;
MEs_new  = MEs_new[,-c(1:2)]

MEs_noGray = MEs_new[,-which(colnames(MEs_new)=="ME0")]

write.csv(MEs_noGray, file="MEs_noGray.csv")
