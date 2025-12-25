if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MAST")

library(MAST)
library(dplyr)
library(Seurat)
library(ggplot2)
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")
BiocManager::install("edgeR")
library(edgeR)


scRNA<- readRDS("chicken_human_.rds")
DefaultAssay(scRNA) <- "RNA"

a <- grep("^RP[SL]",rownames(scRNA))
no.rb.genes <- rownames(scRNA)[-a]
################
scRNA<-subset(x = scRNA, features=no.rb.genes)

data_chi= scRNA[,scRNA@meta.data$species %in% c("chicken")]
data_hum= scRNA[,scRNA@meta.data$species %in% c("human1")]



matrix = as.matrix(GetAssayData(data_hum, slot="counts"))
hum_matrix <- edgeR::cpm(matrix)
matrix = as.matrix(GetAssayData(data_chi, slot="counts"))
chi_matrix <- edgeR::cpm(matrix)

chi_matrix <- log(chi_matrix+1)
hum_matrix <- log(hum_matrix+1)

c <- CreateSeuratObject(count=chi_matrix,meta.data = data_chi@meta.data)
h <- CreateSeuratObject(count=hum_matrix,meta.data = data_hum@meta.data)

merged.anno.comp.ds <- merge(h, y=c)
Idents(merged.anno.comp.ds) <- merged.anno.comp.ds$species

head2head2 <- FindMarkers(merged.anno.comp.ds, ident.1 = "chicken", ident.2 = "human1",test.use = 'MAST',min.pct = 0.25, logfc.threshold = 1)
head2head2$gene <- rownames(head2head2)

head2head2_cpm = head2head2 %>% dplyr::select(gene,everything()) %>% subset(p_val_adj<0.05) %>% subset(abs(avg_log2FC) > 2)

var_featrue = intersect(rownames(data_chi@assays$RNA), rownames(data_hum@assays$RNA))
GetLogUMI <- function(rds_object, gene_list) {
  matrix = as.matrix(GetAssayData(rds_object, slot="counts"))
  matrix = edgeR::cpm(matrix)
  matrix = matrix[gene_list, ]
  matrix = t(matrix)
  x = c()
  for (i in var_featrue){
    cell_num = sum(matrix[, i] > -1)
    gene_sum = sum(matrix[, i])
    value = log(gene_sum/cell_num + 1)
    x = append(x, value)
  }
  return(x)
}

x = GetLogUMI(data_chi, var_featrue)
y = GetLogUMI(data_hum, var_featrue)

data_plot = data.frame(var_featrue, x, y)
colnames(data_plot) = c("gene", "x", "y")
data_plot[is.na(data_plot)] = 0





correlation_matrix <- cor.test(data_plot$x, data_plot$y)
r <- correlation_matrix$estimate
r


max(data_plot$x)
max(data_plot$y)



write.csv(head2head,"chi_human_head2head.csv")


data_plot %>% mutate(size=ifelse(gene %in% head2head2_cpm$gene, "2", "2")) -> data_plot 

data_plot %>% mutate(lable=ifelse(gene %in% head2head2_cpm$gene, "mark", "normal")) -> data_plot 

data_plot2=data_plot[which(data_plot$lable == "mark"), ]

p <- ggplot(data = data_plot,aes(x = x, y = y)) +
  geom_point(aes(),alpha=0.3,size=1.5, color="grey") +
  scale_x_continuous("Chicken") +
  scale_y_continuous("Human")+
  theme_bw()+
  theme(
    legend.background=element_blank(), legend.key=element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  ) +coord_cartesian(xlim = c(0, 12), ylim = c(0,12), expand = F)+
  geom_point(data=data_plot2, mapping=aes(x,y), size=1.5, color="red",alpha=0.5)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = "R=0.57")
p
