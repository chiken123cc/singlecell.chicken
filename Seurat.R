

library(Seurat)


pbmc.data<-Read10X(data.dir ="BM",gene.column = 1) 
pbmc.data<-Read10X(data.dir ="filtered_feature_bc_matrix") 

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "BM_34w_1", min.cells = 3, min.features = 200)

mt.genes <- c('ATP6','ATP8','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5','ND6')


mt.genes2 <- mt.genes[which(mt.genes %in% row.names(pbmc))]

pbmc <- PercentageFeatureSet(object = pbmc, features = mt.genes, col.name = "percent.mt")

pbmcmeta <- pbmc@meta.data[order(-pbmc@meta.data$nFeature_RNA),]
n95 <- as.numeric(as.integer(nrow(pbmcmeta) * 0.05))
n95features <- as.numeric(pbmcmeta[n95,"nFeature_RNA"])
pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < n95features & percent.mt < 15)



pbmc  <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

ElbowPlot(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:15, min.dist = 0.10)
#pbmc <- RunTSNE(pbmc, dims = 1:15, min.dist = 0.10)

DimPlot(pbmc,reduction = 'umap',label=T) + NoLegend()


library(harmony)

pbmc <- RunHarmony(pbmc, "orig.ident")
pbmc <- FindNeighbors(pbmc,dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, reduction = "harmony",dims = 1:15, min.dist = 0.10)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pbmc.markers = pbmc.markers %>% dplyr::select(gene,everything()) %>% subset(p_val_adj<0.05)
top50 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50,file="BF_34w_DEGs.csv")

