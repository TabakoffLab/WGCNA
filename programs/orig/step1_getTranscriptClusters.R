#Set up workspace;
rm(list=ls())
#the last file is the date (aka version) of the WGCNA (year_month_day);
setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218")
library(WGCNA)
options(stringsAsFactors=FALSE)

#####################################################################
# Step1: filter on probesets present in at least 5% of the samples  #
#####################################################################
#load in DABG values
dabg = read.table(file="/data2/saba/ForPhenoGen/HXB.BXH.Brain.fullPS/dabg.fullPS.HXB_BXH.brain.PhenoGen.txt",sep="\t",header=TRUE, row.names=1)
#load in expression data
expr = read.table(file="/data2/saba/ForPhenoGen/HXB.BXH.Brain.fullPS/rma.fullPS.HXB_BXH.brain.PhenoGen.txt",sep="\t",header=TRUE, row.names=1)

table(dim(dabg)==dim(expr))

presentDABG <- dabg[rowSums(dabg<0.0001)>(ncol(expr)*0.05),]
presentExpr <- expr[rowSums(dabg<0.0001)>(ncol(expr)*0.05),]

#176,230 present probesets

table(nrow(presentDABG)==nrow(presentExpr))
table(ncol(presentDABG)==ncol(presentExpr))

save(presentDABG, presentExpr, file="presentDABGandExpr.Rdata")

######################################################
# Step2: find those present probesets within an exon #
######################################################

load("presentDABGandExpr.Rdata")
gtfFile = read.table(file="/data2/saba/BNLx.SHR/RNA-Seq.Brain.total/reconstruction/reconstruct.total.brain.FINAL.26Aug14.gtf",sep="\t", header=FALSE)

gtfFile = gtfFile[,c(1,4,5, 7, 9)]
colnames(gtfFile) = c("chr", "start", "stop", "strand", "name")
gene = sapply(strsplit(sapply(strsplit(as.matrix(gtfFile[,"name"]), split="gene_id ", fixed=TRUE), "[[",2), ";", fixed=TRUE), "[[", 1)
gtfFile = cbind(gtfFile, gene)
gtfFile = gtfFile[,-which(colnames(gtfFile)=="name")]

gtfFile_neg = gtfFile[which(gtfFile[,"strand"]=="-"),]
gtfFile_pos = gtfFile[which(gtfFile[,"strand"]=="+"),]
gtfFile_noStrand = gtfFile[which(gtfFile[,"strand"]=="."),]

###Use BEDTOOLS to find what RNAseq reconstructed transcriptome exons overlap with affymetrix exon probesets
###Need to get data in correct format;

forBEDtools = function(x){
	x2 = x[,-which(colnames(x)=="strand")]
	return(x2)}

gtfFile_neg = forBEDtools(gtfFile_neg)
gtfFile_pos = forBEDtools(gtfFile_pos)
gtfFile_noStrand = forBEDtools(gtfFile_noStrand)

#these are the RNAseq files you will use
write.table(gtfFile_neg, file="gtfFile_neg.txt",sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
write.table(gtfFile_pos, file="gtfFile_pos.txt",sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
write.table(gtfFile_noStrand, file="gtfFile_noStrand.txt",sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)

PSloc = read.table(file="/data2/kiemele/HXB/UpdateRatGenome/try20121012.rn5/results/probeSetLocations.wStrand.ratExonArray.18Apr2013.bed", sep="\t", header=FALSE)
colnames(PSloc) = c("chr", "start", "stop", "probesetID", "strand")

positivePS = PSloc[which(PSloc[,"strand"]=="+"),]
negativePS = PSloc[which(PSloc[,"strand"]=="-"),]

PSall = forBEDtools(PSloc)
PSpos = forBEDtools(positivePS)
PSneg = forBEDtools(negativePS)

#these are the exon array files you will use
write.table(PSall, file="PSall.txt",sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
write.table(PSneg, file="PSneg.txt",sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
write.table(PSpos, file="PSpos.txt",sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)


###Use BEDtools to find overlap between:
# 1. negative PS and negative RNAseq
# 2. positive PS and positive RNAseq
# 3. all PS and no strand identified RNAseq;
# CODE: intersectPSandRNAseqTC.txt

##reload overlapped results into R:
posMatch = read.table(file="positivePSandRNAseqTCoverlap.txt",sep="\t", header=FALSE)
negMatch = read.table(file="negativePSandRNAseqTCoverlap.txt",sep="\t", header=FALSE)
noStrandMatch = read.table(file="noStrandPSandRNAseqTCoverlap.txt",sep="\t", header=FALSE)

getHeaders = function(x){
	x2 = x[,c(4,8)]
	colnames(x2) = c("probeset_id", "gene_id")
	return(x2)
}

posMatchWant = getHeaders(posMatch)
negMatchWant = getHeaders(negMatch)
noStrandMatchWant = getHeaders(noStrandMatch)

exonInfo = rbind(posMatchWant, negMatchWant, noStrandMatchWant)
test = table(exonInfo[,1])
head(test)

length(unique(exonInfo[,"probeset_id"]))
#156,840 unique probesets
length(unique(exonInfo[,"gene_id"]))
#14,695 unique transcript clusters


# The gtf file also had transcript information on it.
# Since we are just using the gene info (so we have no repeated PS in dataset)
# we might have some repeated ps to gene identifiers (if multiple transcripts)
# need to delete any redunant info
exonInfo = exonInfo[order(exonInfo[,"probeset_id"]),]
toDelete = c()
for(i in 2:nrow(exonInfo)){
	toDelete = c(toDelete, exonInfo[i,"probeset_id"]==exonInfo[(i-1),"probeset_id"])
}
length(toDelete)
toDelete = c("FALSE", toDelete)
table(toDelete)

uniqueExonInfo = exonInfo[-which(toDelete=="TRUE"),]

nrow(exonInfo)
length(unique(exonInfo[,"probeset_id"]))
nrow(exonInfo) - length(unique(exonInfo[,"probeset_id"]))

length(unique(uniqueExonInfo[,"probeset_id"]))
dim(uniqueExonInfo)

save(uniqueExonInfo, file="PStoRNAseqGeneID.Rdata")

#intersect the ps located in exons and present ps
presentExonExpr = merge(uniqueExonInfo, presentExpr, by.x="probeset_id", by.y=0)
length(unique(presentExonExpr[,"gene_id"]))
#9,623 unique genes

length(unique(presentExonExpr[,"probeset_id"]))
#82,206 unique probesets
dim(presentExonExpr)

presentKey = merge(uniqueExonInfo, as.matrix(rownames(presentExpr)), by.x="probeset_id", by.y=1)

#############################################
# Step 3: Look for outlier samples (arrays) #
#############################################

sampleTree = flashClust(as.dist(1-cor(presentExonExpr[,-c(1:2)])), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.

pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

##Looks like we have 1 sample outliers (height > 0.03):
## But really just look if there is a clear outlier (I usually don't have a fixed height I use for every single dataset);
toRemove = c("PD_8604_02_Brain.CEL")

sig=c()
for(i in toRemove){
	sig = c(sig, which(colnames(presentExonExpr)==i))
}

length(sig)

PresentExonExpr_clean = presentExonExpr[,-sig] 
rownames(PresentExonExpr_clean) = PresentExonExpr_clean[,1]
PresentExonExpr_clean = PresentExonExpr_clean[,-c(1:2)]

####################################################################################
# Step 4: Correlate strain means & determine what PS should be collapsed together  #
####################################################################################
strain = sapply(strsplit(colnames(PresentExonExpr_clean), split="_", fixed=TRUE), "[[", 1)
means = t(apply(PresentExonExpr_clean, 1, function(a) lm(a~as.factor(strain) -1)$coefficients))
colnames(means) = sapply(strsplit(colnames(means), split=")", fixed=TRUE), "[[", 2)

#remove parents and the other inbreds (PD, WKY.Lj) from correlation analysis
toDelete = c("BN.LX", "PD", "SHR.H", "SHR.lj", "SHR.Lx", "WKY.Lj")
sig = c()
for(i in toDelete){
	sig = c(sig, which(colnames(means)==i))
}

length(sig)
RImeans = means[,-sig]

#make a list of matrices, where each matrix represents the expression matrix for a specific TC
TCmatrices = function(x, key){
	xList = list()
	sig = list()
	TC = as.matrix(unique(key[,"gene_id"]))
	for(i in TC){
		sig[[i]] = key[which(key[,"gene_id"]==i),"probeset_id"]
		xList[[i]] = merge(x, as.matrix(sig[[i]]), by.x=0, by.y=1)
		rownames(xList[[i]]) = xList[[i]]$Row.names
		xList[[i]] = xList[[i]][,-1]
	}
	return(xList)
}

TCmeansMatrices = TCmatrices(RImeans, presentKey)

#now we want to know what ps can be collapsed down together (those with correlation > 0.25)
#then we can make a function and get the 1st principal component when for those with the individual array data;

cluster = function(x){
  if(nrow(x)>1){	
  tree = hclust(as.dist(1-cor(t(x))), method = "average")
  cutTree = cutree(tree, h=0.75)
   }
 if(nrow(x)==1){
  cutTree = 1
  names(cutTree) = rownames(x)
  }
 return(cutTree)		
}


##example of clustering...
sampleTree = hclust(as.dist(1-cor(t(TCmeansMatrices[[1]]))), method = "average")
sizeGrWindow(12,9)
pdf(file = "exampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Clustering Example: gene XLOC_022227", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 0.75, col="red")
abline(h = 0.5, col="blue")

dev.off()

#5 if height = 0.75
#9 if height = 0.5

cluster(TCmeansMatrices[[1]])

##Get a list of genes. Each gene contains all probesets that match and then have a # assigned designating to specific cluster
##This takes awhile to run (maybe an hour or 2);
clustersToMake = lapply(TCmeansMatrices, cluster)

getClustNames = function(x){
	names = list()
	PS = list()
	for(i in 1:length(x)){
	  names[[i]] = as.matrix(paste(names(x)[[i]], ".clust", x[[i]], sep=""))
 	  PS[[i]] = names(x[[i]])	
	}
	collapsedResults = cbind(unlist(PS), unlist(names)) 
	colnames(collapsedResults) = c("probeset", "transcriptClusterCollapsed")
	return(collapsedResults)
}

PSclustKey = getClustNames(clustersToMake)
dim(PSclustKey)

length(unique(PSclustKey[,"transcriptClusterCollapsed"]))
save(PSclustKey, clustersToMake, TCmeansMatrices, file="findingCollapsedTC.Rdata")
write.csv(PSclustKey, file="PStoCollapsedTransciptClusterKey.csv")

################################################################################################
# Step 5: Generate the 1st principal component scores for those that need to be collapsed down #
################################################################################################
#want to use the PresentExonExpr_clean... But need to remove the parentals;

strain = sapply(strsplit(colnames(PresentExonExpr_clean), split="_", fixed=TRUE), "[[", 1)
#remove parents and the other inbreds (PD, WKY.Lj) when generating the PCA analyses;
toDelete = c("BN.LX", "PD", "SHR.H", "SHR.lj", "SHR.Lx", "WKY.Lj")
sig = c()
for(i in toDelete){
	sig = c(sig, which(strain==i))
}

length(sig)
RIforCollapse = PresentExonExpr_clean[,-sig]

TCcollapsed = collapseRows(RIforCollapse, PSclustKey[,"transcriptClusterCollapsed"], PSclustKey[,"probeset"], method="ME", thresholdCombine=NA)
length(TCcollapsed)
save(TCcollapsed, file="TCcollapsed.datExpr.Rdata")
#17,579 TC

##################################################################################################
# Step 6: Calculate the heritability and apply a filter of > 0.25 to stay in dataset for network #
##################################################################################################
expr = TCcollapsed$datETcollapsed
strain = sapply(strsplit(colnames(expr), split="_", fixed=TRUE), "[[", 1)

test = apply(expr, 1, function(a) summary(lm(a~as.factor(strain)))$r.squared)

herits = cbind(rownames(expr), test)
colnames(herits) = c("transcript", "heritability")
write.table(herits, file="herits.clustersTry20141218.HXB.brain.txt", sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

heritable <- herits[herits[,2]>0.25,1]
heritablePS <- match(heritable,rownames(expr))
heritablePS <- heritablePS[!is.na(heritablePS)]
heritExprs <- expr[heritablePS,]
save(heritExprs, file="heritExprs.Rdata")
dim(heritExprs)
#12,450 heritable TC

##look heritability across all transcript clusters
pdf(file = "heritabilities.pdf", width = 12, height = 9);
hist(as.numeric(herits[,2]), breaks=100, main = "Histogram of Heritability", xlim=c(0,1), axes=TRUE)
abline(v = 0.25, col="red")
dev.off()

#get transcript cluster means for the heritable transcript cluster (to use in WGCNA analysis)
heritMeans = t(apply(heritExprs, 1, function(a) lm(a~as.factor(strain) -1)$coefficients))
colnames(heritMeans) = sapply(strsplit(colnames(heritMeans), split=")", fixed=TRUE), "[[", 2)
save(heritMeans, file="heritExprs_strainMeans.Rdata")

####################################
# Step 7: Look for outlier strains #
####################################
sampleTree = flashClust(as.dist(1-cor(heritMeans)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
pdf(file = "strainClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Strain clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
#keep all strains;


