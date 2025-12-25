


#######################################################################tcell.Pseudotime

library(monocle)
library(dplyr)
library(Seurat)

pbmc <- readRDS("pbmc.rds")

expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
p_data <- pbmc@meta.data 
f_data <- data.frame(gene_short_name = row.names(pbmc),
                     row.names = row.names(pbmc))

pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) 
length(expressed_genes)


Time_diff <- differentialGeneTest(cds, cores = 10, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

write.csv(Time_diff,'Time_diff.csv')

gene<-c(
"TCF12", "GATA2", "FOS", "KLF7", "RUNX3", "NME2", "HOXB6", "TCF3", "CREB3L2", "NFATC1", "FOSL2", "CEBPD", "SMAD3", "BATF", "EGR1", "ZNF366", "ERG", "HOXA5", "CREM", "HOXA9",
 "HOXA3", "LITAF", "CEBPB", "NCOR2", "ST18", "JUN", "MAFB", "PLEK", "SPIC", "BCL11A", "MYB")
 
pdf("heatmap.new3.pdf",width =5,height =6)
plot_pseudotime_heatmap(cds[gene,],
                        num_clusters =2,show_rownames =T,
                        cores = 15,hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
dev.off()


plot_pseudotime_heatmap(cds[gene,],
                        ##add_annotation_col = your_column,
                        cluster_rows = TRUE,
                        trend_formula = "~sm.ns(Pseudotime, df=1)",
                        hclust_method = ALL_HC, 
                        num_clusters = 1,
                        hmcols = my_color_set[sub_color][[1]],
                        scale_max = 3, 
                        scale_min = -3,
                        cores = 20,
                        show_rownames = T,
                        return_heatmap = FALSE)
												
