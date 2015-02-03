
#################################################################################################
# Get adjacency matrices to find out connectivity for all TC within candiate modules & hub gene #
##################################################################################################
setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218")
library(WGCNA)
options(strainsAsFactors=FALSE)

load("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141218/heritExprs_strainMeans.Rdata")

###############################
# Create a Correlation Matrix #
###############################

corMatrix = bicor(t(heritMeans))
dim(corMatrix)

###############################
# Create the Adjacency Matrix #
###############################

adjMatrix_all = corMatrix^9
dim(adjMatrix_all)

save(adjMatrix_all, file="adjMatrix_all.Rdata")

###################################################
# Create a List of Adjacency Matrices (by Module) #
###################################################
key = read.csv("moduleMembership.csv")

table(rownames(adjMatrix_all)==key[,"TC"])

modLabs = names(table(key[,"moduleLabels"]))
modLabs = modLabs[-which(modLabs=="0")]

TCwant = list()
adjMatricesByModule = list()
for(i in modLabs){	
	TCwant[[i]] = as.vector(key[which(key[,"moduleLabels"]==i),"TC"])
	adjMatricesByModule[[i]] = adjMatrix_all[TCwant[[i]], TCwant[[i]]]
}

names(adjMatricesByModule)
key2 = unique(key[,c(2:3)])

test = c()
for(i in names(adjMatricesByModule)){
	test = c(test, as.character(key2[which(key2[,"moduleLabels"]==i), "moduleColors"]))
}

names(adjMatricesByModule) = test
save(adjMatricesByModule, file="adjMatricesByModule.Rdata")


