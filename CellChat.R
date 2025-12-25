library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
data.input <- pbmc@assays$RNA@data
meta <- pbmc@meta.data
meta$labels <- pbmc$cell_type
cell.use = rownames(meta)[meta$symbol == "young"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels")
CellChatDB <- CellChatDB.human（mouse）
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
