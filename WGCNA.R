
################################################


library(WGCNA)

pbmc<-readRDS("module.t.tur.chi.hum.rds")

pbmc@meta.data[["cell_type"]]<-Idents(pbmc)

###########1
datadf <- as.matrix(GetAssayData(pbmc, slot = "data"))
idd1 <- pbmc@meta.data
Inter.id1 <- data.frame(Inter.id = rownames(idd1), cell_type = as.character(idd1$cell_type))
rownames(Inter.id1)<-rownames(idd1)
colnames(Inter.id1)<-c("CellID","Celltype")
Inter.id1<-as.data.frame(Inter.id1)
head(Inter.id1)
Inter1<-datadf[,Inter.id1$CellID]
Inter2<-as.matrix(Inter1)
Inter2[1:4,1:4]


##########
new_ids_list1 = list()


for (i in 1:length(levels(pbmc))) {
  cluster_id = levels(pbmc)[i]
  cluster_cells <- rownames(Inter.id1[Inter.id1$Celltype == cluster_id,])
  cluster_size <- length(cluster_cells)     
  pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
  pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
  names(pseudo_ids) <- sample(cluster_cells)    
  new_ids_list1[[i]] <- pseudo_ids      
}



new_ids <- unlist(new_ids_list1)
new_ids <- as.data.frame(new_ids)
new_ids_length <- table(new_ids)


new_colnames <- rownames(new_ids)  
all.data<-datadf[,as.character(new_colnames)] 
all.data <- t(all.data)
new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                    list(name=new_ids[,1]),FUN=mean)
rownames(new.data)<-new.data$name
new.data<-new.data[,-1]
new_ids_length<-as.matrix(new_ids_length)
head(new.data)[1:5,1:5]

#######################################		
#softPower
dataExpr <- as.data.frame(new.data)			   
SubGeneNames=colnames(dataExpr)

powers = c(c(1:18), seq(from = 20, to=40, by=2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5,networkType = "signed")

pdf("p1.pdf", width=12)
par(mfrow = c(1,2));
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()




softPower = 10
NetworkType = "signed"
pbmc_adj= adjacency(dataExpr,type = NetworkType, power = softPower)
pbmc_k <- softConnectivity(dataExpr, type="signed", power=softPower)
pbmc_k <- pbmc_k/max(pbmc_k)
names(pbmc_k) <- colnames(dataExpr)
quantile(pbmc_k, probs=seq(0,1,by=0.1))




######## find modules pbmc ######
###################################
library(flashClust)	

pbmc_TOM=TOMsimilarityFromExpr(dataExpr,networkType = NetworkType, TOMType = "signed", power = softPower)
colnames(pbmc_TOM) =rownames(pbmc_TOM) =SubGeneNames

pbmc_dissTOM = 1-pbmc_TOM
pbmc_geneTree = flashClust(as.dist(pbmc_dissTOM),method="average")

pdf("p2.pdf", width=12)
plot(pbmc_geneTree, xlab="", sub="",cex=0.3, labels=F)
dev.off()

pbmc_minModuleSize = 30
pbmc_dynamicMods = cutreeDynamic(dendro = pbmc_geneTree, distM = pbmc_dissTOM, method="hybrid", deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = pbmc_minModuleSize, pamStage=T)
pbmc_dynamicColors = labels2colors(pbmc_dynamicMods)
table(pbmc_dynamicColors)
pbmc_dynamicColors

  

plotDendroAndColors(pbmc_geneTree, pbmc_dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")


# diag(pbmc_dissTOM) = NA
# module merging

pbmc_MEList = moduleEigengenes(dataExpr, colors = pbmc_dynamicColors)
pbmc_MEs = pbmc_MEList$eigengenes

pbmc_MEDiss = 1-cor(pbmc_MEs)
pbmc_METree = hclust(as.dist(pbmc_MEDiss), method = "average")

plot(pbmc_METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")



pbmc_MEDissThres = 0.5
# Call an automatic merging function
pbmc_merge = mergeCloseModules(dataExpr, pbmc_dynamicColors, cutHeight = pbmc_MEDissThres, verbose = 3)
# The merged module colors
pbmc_mergedColors = pbmc_merge$colors
# Eigengenes of the new merged modules:
pbmc_mergedMEs = pbmc_merge$newMEs


plotDendroAndColors(pbmc_geneTree, cbind(pbmc_dynamicColors, pbmc_mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)




pbmc_dynamicColors <- pbmc_mergedColors

##########
library(dplyr)

pbmc_MEList = moduleEigengenes(dataExpr, colors = pbmc_dynamicColors)
pbmc_MEs = pbmc_MEList$eigengenes

pbmc_MEeigengenes <- as.matrix(pbmc_MEList$eigengenes)
rownames(pbmc_MEeigengenes) <- rownames(dataExpr)



pbmc_tempident <- rownames(pbmc_MEeigengenes)
pbmc_tempident <- sapply(pbmc_tempident, function(x) gsub("_Cell[0-9]{1,2}$", "", x))
pbmc_tempident

names(pbmc_tempident) <- rownames(pbmc_MEeigengenes)
pbmc_tempident <- as.factor(pbmc_tempident)



pbmc_tempident <- as.data.frame(pbmc_tempident)
pbmc_tempident <- cbind(pbmc_tempident, rownames(pbmc_tempident))
colnames(pbmc_tempident) <- c("tempident","pbmc_id")
pbmc_tempident <- arrange(pbmc_tempident, tempident)
pbmc_tempident2 <- as.factor(pbmc_tempident[,1])
names(pbmc_tempident2) <- pbmc_tempident[,2]
pbmc_tempident <- pbmc_tempident2

library(RColorBrewer)
library(plyr)

n.clusters = length(table(pbmc_tempident))
plot_colors = colorRampPalette(brewer.pal(9, "Set1"))(n.clusters)
pbmc_plotColors <- mapvalues(pbmc_tempident, from = levels(pbmc_tempident), to = plot_colors) # plyr package
pbmc_plotColors <- as.character(pbmc_plotColors)

###########################
pbmc_MEeigengenes <- as.matrix(pbmc_MEList$eigengenes)
pbmc_MEeigengenes[pbmc_MEeigengenes>0.3] <- 0.3
rownames(pbmc_MEeigengenes) <- rownames(dataExpr)
pbmc_MEeigengenes <- pbmc_MEeigengenes[order(match(rownames(pbmc_MEeigengenes), names(pbmc_tempident))),]


library('gplots')

par(oma=c(4,4,4,10))
colors <- colorRampPalette(c("#9D8BFF","black", "#FFE000"))(n = 200)
pal <- colorRampPalette(colors)
heatmap.2(t(pbmc_MEeigengenes), Colv=F, Rowv=F, trace="none",  col=pal, ColSideColors= pbmc_plotColors,na.color="green")
dev.off()



# extract names of genes that belong to the modules
SubGeneNames=colnames(dataExpr)
pbmc_module_colors= setdiff(unique(pbmc_dynamicColors), "grey")
for (color in pbmc_module_colors){
    module=SubGeneNames[which(pbmc_dynamicColors==color)]
    write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)    
}


	