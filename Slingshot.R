
library(tradeSeq)
library(slingshot)
library(SingleCellExperiment)
library(qs)
library(tidyverse)
library(RColorBrewer)



scobj = readRDS(file = 'pbmc2.rds')
counts <- scobj@assays$RNA@counts
sim <- SingleCellExperiment(assays = List(counts = counts)) 
umap = scobj@reductions$umap@cell.embeddings
colnames(umap) = c('UMAP-1', 'UMAP-2')
reducedDims(sim) = SimpleList(UMAP = umap)
meta = scobj@meta.data
colData(sim)$cell_type = meta$cell_type


sim <- slingshot(sim, 
                 clusterLabels = 'cell_type',  
                 reducedDim = 'UMAP',  
                 start.clus= "tn",  
                 end.clus = "tx"   
)
   
colnames(colData(sim))



colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) 
plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)] 
plotcol[is.na(plotcol)] <- "lightgrey"



plot(reducedDims(sim)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))


#############


counts <- sim@assays@data$counts
crv <- SlingshotDataSet(sim)

set.seed(111)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
system.time({ 
  sce <- fitGAM(counts = counts, 
                pseudotime = pseudotime, 
                cellWeights = cellWeights,
                nknots = 6, 
                verbose = FALSE)
})
table(rowData(sce)$tradeSeq$converged)


topgenes <- marker$x
			  		  
			  
pst.ord <- order(sim$slingPseudotime_1, na.last = NA)
heatdata <- assays(sim)$counts[topgenes, pst.ord]
heatclus <- sim$GMM[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

