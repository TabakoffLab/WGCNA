######################################################
##  Creating New Transcript Clusters: ALL MPS File  ##
##  Version 11/13/2014	     					##
##  Note: This file will be for finding the DABG in ##
##  the RNAseq derived TC.  This is for the WGCNA   ##
##  results posted on Phenogen ONLY.  Users can see ##
##  easily what is above and below background, but  ##
##  this is NOT for use in filtering prior to TC    ##
##  construction. 									##
######################################################


### Load in PS to TC key:
setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/dabgForPhenogen")
rm(list=ls())
##Load in Final Set of Transcript Clusters used for analysis
load(file="/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/heritExprs_strainMeans.Rdata")
##Load in File with the PS to TC key;
load(file="/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/findingCollapsedTC.Rdata")

PStoTCkey_usedInAnalysis = merge(as.matrix(rownames(heritMeans)), PSclustKey, by.x=1, by.y="transcriptClusterCollapsed")
colnames(PStoTCkey_usedInAnalysis) = c("transcript_cluster_id", "probeset_id")

dim(PStoTCkey_usedInAnalysis)
#62,103 probesets
length(unique(PStoTCkey_usedInAnalysis[,"transcript_cluster_id"]))
#12,450 unqiue TC


setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/dabgForPhenogen")
### Load in file With PS and # probes (Spencer generated this file);
psToNumProbes = read.table(file="HXB.rn5.probe_count.txt", sep="\t", header=TRUE)
colnames(psToNumProbes) = c("probeset_id", "probe_count")

### Merge these files into 1;
RNAseqTC.mps = merge(PStoTCkey_usedInAnalysis, psToNumProbes, by="probeset_id")
dim(RNAseqTC.mps)
RNAseqTC.mps = RNAseqTC.mps[order(RNAseqTC.mps[,"transcript_cluster_id"]),]
write.table(RNAseqTC.mps, file="RNAseqTC.mps", sep="\t", row.names=FALSE, quote=FALSE)


