library(future)
plan("multisession", workers = 12)
options(future.globals.maxSize= 62914560000)
library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)



one2one_ch <- read.csv("one2one_gene.csv")


chi <- readRDS("Chicken.rds")
chi_raw_data <- chi@assays$RNA@counts
chi_raw_data <- chi_raw_data[which(row.names(chi_raw_data) %in%  one2one_ch$chicken_gene),]
row.names(chi_raw_data) <- one2one_ch$human_gene[match(row.names(chi_raw_data), one2one_ch$chicken_gene)]


hum <- readRDS("Human.rds")
hum_raw_data <- hum@assays$RNA@counts



common_genes <- intersect(row.names(chi_raw_data),row.names(hum_raw_data))
chi_raw_data <- chi_raw_data[match (common_genes,row.names(chi_raw_data)),]
hum_raw_data <- hum_raw_data[match (common_genes,row.names(hum_raw_data)),]


chi_new <- CreateSeuratObject(chi_raw_data,project = "stoma",meta.data = chi@meta.data)
hum_new <- CreateSeuratObject(hum_raw_data,project = "stoma",meta.data = hum@meta.data)
chi_new$species <- "Chicken"
hum_new$species <- "Human"
pbmc.list <- list(chi_new,hum_new)

saveRDS(pbmc.list,"stoma_chi_hum.rds")


for(i in 1:length(pbmc.list)) {
  pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], return.only.var.genes = FALSE)
}

nVarGenes <- 1500
nDims_CCA <- 40

int_features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = nVarGenes)
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = int_features,  verbose = FALSE)
pbmc_anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features= int_features, normalization.method = "SCT", dims = 1:40,reduction="cca")
pbmc <- IntegrateData(anchorset = pbmc_anchors, normalization.method = "SCT", dims = 1:40)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
dim.usage <- 18
pbmc <- FindNeighbors(pbmc, dims = 1:dim.usage)
pbmc <- FindClusters(pbmc, resolution = 0.6)
pbmc <- RunUMAP(pbmc, dims = 1:dim.usage, min.dist = 0.3)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc,file="chi_hum.rds")


