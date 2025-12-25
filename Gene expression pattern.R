
##############################################################gene.module.evo



library(Seurat)
library(tidyverse)
library(viridis)
library(RColorBrewer)


htiss.anno.comp<-readRDS('human.rds')

mtiss.anno.comp<-readRDS('chicken.rds')

#################

celltypes<-levels(htiss.anno.comp)
ngenes <- length(rownames(htiss.anno.comp@assays$RNA@data))
cnum <- length(celltypes)

average_cells.comp <- matrix(0,ngenes,cnum)
pct_exp.comp <- matrix(0,ngenes,cnum)

colnames(average_cells.comp) <- celltypes
rownames(average_cells.comp) <- rownames(htiss.anno.comp@assays$RNA@data)

colnames(pct_exp.comp) <- celltypes
rownames(pct_exp.comp) <- rownames(htiss.anno.comp@assays$RNA@data)


for (i in 1:cnum) {
  cells <- WhichCells(htiss.anno.comp, ident = celltypes[i])
  average_cells.comp[,i] <- rowMeans(as.matrix(htiss.anno.comp@assays$RNA@data[,cells]))
  average_cells.comp[which(is.nan(average_cells.comp))] <- 0
  pct_exp.comp[,i] <- rowSums(htiss.anno.comp@assays$RNA@data[,cells] > 0) / length(cells)
}


jitter <- c(-2, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2)
is_on.comp <- list()


for (j in jitter) {
  for (k in jitter) {
    
    tmp <- matrix(0,ngenes,cnum)
    
    colnames(tmp) <- celltypes
    rownames(tmp) <- rownames(htiss.anno.comp@assays$RNA@data)
    
    for (i in 1:cnum) {
      
      median_cutoff <- median(average_cells.comp[,i]) + j * sd(average_cells.comp[,i])
      pct_cutoff <- median(pct_exp.comp[,i]) + k * sd(pct_exp.comp[,i])
      
      tmp[intersect(names(which(average_cells.comp[,i] > median_cutoff)), names(which(pct_exp.comp[,i] > pct_cutoff))),i] <- 1
      
    }
    
    colnames(tmp) <- celltypes
    rownames(tmp) <- rownames(htiss.anno.comp@assays$RNA@data)
    
    is_on.comp[[paste0(j,'-',k)]] <- tmp
    
  }
}


#################################

chicken.celltypes<-levels(mtiss.anno.comp)
ngenes <- length(rownames(mtiss.anno.comp@assays$RNA@data))
chicken.cnum<-length(chicken.celltypes)





chicken.average_cells.comp <- matrix(0,ngenes,chicken.cnum)
chicken.pct_exp.comp <- matrix(0,ngenes,chicken.cnum)

colnames(chicken.average_cells.comp) <- chicken.celltypes
rownames(chicken.average_cells.comp) <- rownames(mtiss.anno.comp@assays$RNA@data)

colnames(chicken.pct_exp.comp) <- chicken.celltypes
rownames(chicken.pct_exp.comp) <- rownames(mtiss.anno.comp@assays$RNA@data)

for (i in 1:chicken.cnum) {
  cells <- WhichCells(mtiss.anno.comp, ident = chicken.celltypes[i])
  chicken.average_cells.comp[,i] <- rowMeans(as.matrix(mtiss.anno.comp@assays$RNA@data[,cells]))
  chicken.average_cells.comp[which(is.nan(chicken.average_cells.comp))] <- 0
  chicken.pct_exp.comp[,i] <- rowSums(mtiss.anno.comp@assays$RNA@data[,cells] > 0) / length(cells)
  
}


##########################
chicken.is_on.comp <- list()

for (j in jitter) {
  for (k in jitter) {
    
    tmp <- matrix(0,ngenes,chicken.cnum)
    
    colnames(tmp) <- chicken.celltypes
    rownames(tmp) <- rownames(mtiss.anno.comp@assays$RNA@data)
    
    for (i in 1:chicken.cnum) {
      
      median_cutoff <- median(chicken.average_cells.comp[,i]) + j * sd(chicken.average_cells.comp[,i])
      pct_cutoff <- median(chicken.pct_exp.comp[,i]) + k * sd(chicken.pct_exp.comp[,i])
      
      tmp[intersect(names(which(chicken.average_cells.comp[,i] > median_cutoff)), names(which(chicken.pct_exp.comp[,i] > pct_cutoff))),i] <- 1
      
    }
    
    colnames(tmp) <- chicken.celltypes
    rownames(tmp) <- rownames(mtiss.anno.comp@assays$RNA@data)
    
    chicken.is_on.comp[[paste0(j,'-',k)]] <- tmp
    
  }
}





class_numbers_list <- list()
class_1_chicken_list <- list()
class_1_human_list <- list()
class_2_chicken_list <- list()
class_2_human_list <- list()
class_2_both_list <- list()
class_3_list <- list()
conserved_exp_list <- list()
conserved_not_exp_list <- list()


chicken.order <- c(1:22)
human.order <- c(1:22)


for (j in jitter) {
  print(j)
  for (k in jitter) {
    
    class_numbers <- matrix(0,ngenes,12)
    
    rownames(class_numbers) <- rownames(htiss.anno.comp@assays$RNA@counts)
    
    class_1_chicken <- c()
    class_1_human <- c()
    class_2_chicken <- c()
    class_2_human <- c()
    class_2_both <- c()
    class_3 <- c()
    conserved_exp <- c()
    conserved_not_exp <- c()
    
    for (i in rownames(htiss.anno.comp@assays$RNA@data)) {
      
      a <- is_on.comp[[paste0(j,'-',k)]][i,][human.order]
      b <- chicken.is_on.comp[[paste0(j,'-',k)]][i,][chicken.order]
      
      ACC <- (a - b) + (2 * (a * b))
      
      # a    1  0  1  1  0
      # b    1  0  0  1  1
      # ACC  2  0  1  2 -1
      
      
      # Class 2 - Expansion & Conserved Exp
      if (any(ACC == 2)) {
        
        if (any(ACC == 1) && !any(ACC == -1)) {
          class_2_human <- c(class_2_human, i)
          class_numbers[i,3] <- length(which(ACC == 2))
          class_numbers[i,4] <- length(which(ACC == 1))
          
        } else if (!any(ACC == 1) && any(ACC == -1)) {
          class_2_chicken <- c(class_2_chicken, i)
          class_numbers[i,5] <- length(which(ACC == 2))
          class_numbers[i,6] <- length(which(ACC == -1))
          
        } else if (any(ACC == 1) && any(ACC == -1)) {
          class_2_both <- c(class_2_both, i)
          class_numbers[i,7] <- length(which(ACC == 2))
          class_numbers[i,8] <- length(which(ACC == 1))
          class_numbers[i,9] <- length(which(ACC == -1))
          
        } else {
          conserved_exp <- c(conserved_exp, i)
          class_numbers[i,12] <- length(which(ACC == 2))
          
        }
        
        # Not Expressed
      } else if (all(ACC == 0)) {
        conserved_not_exp <- c(conserved_not_exp, i)
        
        # Class 3 - Switching
      } else if (any(sign(ACC) > 0) && any(sign(ACC) < 0)) {
        class_3 <- c(class_3, i)
        class_numbers[i,10] <- length(which(ACC == 1))
        class_numbers[i,11] <- length(which(ACC == -1))
        
        # Class 1 - Loss / Gain
      } else if (any(sign(ACC) > 0) && !any(sign(ACC) < 0)) {
        class_1_human <- c(class_1_human,i)
        class_numbers[i,1] <- length(which(ACC == 1))
        
      } else if (!any(sign(ACC) > 0) && any(sign(ACC) < 0)) {
        class_1_chicken <- c(class_1_chicken,i)
        class_numbers[i,2] <- length(which(ACC == -1))
      }
      
    }
    
    class_numbers_list[[paste0(j,'-',k)]] <- class_numbers
    class_1_chicken_list[[paste0(j,'-',k)]] <- class_1_chicken
    class_1_human_list[[paste0(j,'-',k)]] <- class_1_human
    class_2_chicken_list[[paste0(j,'-',k)]] <- class_2_chicken
    class_2_human_list[[paste0(j,'-',k)]] <- class_2_human
    class_2_both_list[[paste0(j,'-',k)]] <- class_2_both
    class_3_list[[paste0(j,'-',k)]] <- class_3
    conserved_exp_list[[paste0(j,'-',k)]] <- conserved_exp
    conserved_not_exp_list[[paste0(j,'-',k)]] <- conserved_not_exp
    
  }
}




#########################################



conserved_exp_avg <- matrix(0,17,17)
colnames(conserved_exp_avg) <- jitter
rownames(conserved_exp_avg) <- jitter

for (j in jitter) {
  for (k in jitter) {
    tmp <- conserved_exp_list[[paste0(j,'-',k)]]
    conserved_exp_avg[as.character(j),as.character(k)] <- length(tmp)
  }
}


library(pheatmap)
pheatmap(conserved_exp_avg)


##
conserved_not_exp_avg <- matrix(0,17,17)
colnames(conserved_not_exp_avg) <- jitter
rownames(conserved_not_exp_avg) <- jitter

for (j in jitter) {
  for (k in jitter) {
    tmp <- conserved_not_exp_list[[paste0(j,'-',k)]]
    conserved_not_exp_avg[as.character(j),as.character(k)] <- length(tmp)
  }
}



class_1_human_avg <- matrix(0,17,17)
colnames(class_1_human_avg) <- jitter
rownames(class_1_human_avg) <- jitter

for (j in jitter) {
  for (k in jitter) {
    tmp <- class_1_human_list[[paste0(j,'-',k)]]
    class_1_human_avg[as.character(j),as.character(k)] <- length(tmp)
  }
}


class_1_chicken_avg <- matrix(0,17,17)
colnames(class_1_chicken_avg) <- jitter
rownames(class_1_chicken_avg) <- jitter

for (j in jitter) {
  for (k in jitter) {
    tmp <- class_1_chicken_list[[paste0(j,'-',k)]]
    class_1_chicken_avg[as.character(j),as.character(k)] <- length(tmp)
  }
}



class_2_human_avg <- matrix(0,17,17)
colnames(class_2_human_avg) <- jitter
rownames(class_2_human_avg) <- jitter

for (j in jitter) {
  for (k in jitter) {
    tmp <- class_2_human_list[[paste0(j,'-',k)]]
    class_2_human_avg[as.character(j),as.character(k)] <- length(tmp)
  }
}


class_2_chicken_avg <- matrix(0,17,17)
colnames(class_2_chicken_avg) <- jitter
rownames(class_2_chicken_avg) <- jitter

for (j in jitter) {
  for (k in jitter) {
    tmp <- class_2_chicken_list[[paste0(j,'-',k)]]
    class_2_chicken_avg[as.character(j),as.character(k)] <- length(tmp)
  }
}



class_2_both_avg <- matrix(0,17,17)
colnames(class_2_both_avg) <- jitter
rownames(class_2_both_avg) <- jitter

for (j in jitter) {
  for (k in jitter) {
    tmp <- class_2_both_list[[paste0(j,'-',k)]]
    class_2_both_avg[as.character(j),as.character(k)] <- length(tmp)
  }
}


class_3_avg <- matrix(0,17,17)
colnames(class_3_avg) <- jitter
rownames(class_3_avg) <- jitter

for (j in jitter) {
  for (k in jitter) {
    tmp <- class_3_list[[paste0(j,'-',k)]]
    class_3_avg[as.character(j),as.character(k)] <- length(tmp)
  }
}

library(pheatmap)
pheatmap(class_3_avg,scale = "none",cluster_rows = F,cluster_cols = F)
pheatmap(conserved_exp_avg,scale = "none",cluster_rows = F,cluster_cols = F)
pheatmap(conserved_not_exp_avg,scale = "none",cluster_rows = F,cluster_cols = F)
pheatmap(class_3_avg,scale = "none",cluster_rows = F,cluster_cols = F)



library(gplots)
heatmap.2(as.matrix(conserved_exp_avg[1:9,9:17]),
          symbreaks = F, symkey = F,  
          key.title = '',
          col=colorRampPalette(brewer.pal(9,'Greys')[2:9])(100),
          trace="none", Rowv = F, Colv = F,
          key.xlab = 'Genes', key.ylab = '',
          density.info = 'density',
          denscol = 'black', labCol = jitter[9:17], 
          labRow = jitter[1:9],
          margins = c(10,10), cexRow = 1)

heatmap.2(as.matrix(conserved_not_exp_avg[1:9,9:17]), symbreaks = F, symkey = F,  key.title = '', col=colorRampPalette(brewer.pal(9,'Greys')[2:9])(100), trace="none", Rowv = F, Colv = F, key.xlab = 'Genes', key.ylab = '', density.info = 'density', denscol = 'black', labCol = jitter[9:17], labRow = jitter[1:9], margins = c(10,10), cexRow = 1)
heatmap.2(as.matrix(class_1_human_avg[1:9,9:17]), symbreaks = F, symkey = F,  key.title = '', col=colorRampPalette(brewer.pal(9,'Greys')[2:9])(100), trace="none", Rowv = F, Colv = F, key.xlab = 'Genes', key.ylab = '', density.info = 'density', denscol = 'black', labCol = jitter[9:17], labRow = jitter[1:9], margins = c(10,10), cexRow = 1)

median(conserved_not_exp_avg)
mean(conserved_not_exp_avg)

install.packages('plotrix')
library('plotrix')
vector_data <- c(conserved_not_exp_avg)
std.error(vector_data)

median(conserved_exp_avg)
mean(conserved_exp_avg)

install.packages('plotrix')
library('plotrix')
vector_data <- c(conserved_exp_avg)
std.error(vector_data)



median(class_3_avg)
mean(class_3_avg)


median(class_1_human_avg)
median(class_1_chicken_avg)


save.image('gene.evo.rdata')
#####

Table_S7 <- matrix(0,ngenes,5)
rownames(Table_S7) <- rownames(htiss.anno.comp@assays$RNA@data)
colnames(Table_S7) <- c('Evo_Type', 'Gene_Class', 'Conserved_Num', 'Human_Spec_Num', 'chicken_Spec_Num')

Table_S7[conserved_exp_list[["0-1"]],'Evo_Type'] <- 0
Table_S7[conserved_not_exp_list[["0-1"]],'Evo_Type'] <- NA
Table_S7[class_1_human_list[["0-1"]],'Evo_Type'] <- 1
Table_S7[class_1_chicken_list[["0-1"]],'Evo_Type'] <- 1
Table_S7[class_2_human_list[["0-1"]],'Evo_Type'] <- 2
Table_S7[class_2_chicken_list[["0-1"]],'Evo_Type'] <- 2
Table_S7[class_2_both_list[["0-1"]],'Evo_Type'] <- 2
Table_S7[class_3_list[["0-1"]],'Evo_Type'] <- 3




Table_S7[names(class_numbers_list[["0-1"]][class_1_human_list[["0-1"]],1]),'Human_Spec_Num'] <- class_numbers_list[["0-1"]][class_1_human_list[["0-1"]],1]
Table_S7[names(class_numbers_list[["0-1"]][class_1_chicken_list[["0-1"]],2]),'chicken_Spec_Num'] <- class_numbers_list[["0-1"]][class_1_chicken_list[["0-1"]],2]

Table_S7[names(class_numbers_list[["0-1"]][class_2_human_list[["0-1"]],3]),'Conserved_Num'] <- class_numbers_list[["0-1"]][class_2_human_list[["0-1"]],3]
Table_S7[names(class_numbers_list[["0-1"]][class_2_human_list[["0-1"]],4]),'Human_Spec_Num'] <- class_numbers_list[["0-1"]][class_2_human_list[["0-1"]],4]

Table_S7[names(class_numbers_list[["0-1"]][class_2_chicken_list[["0-1"]],5]),'Conserved_Num'] <- class_numbers_list[["0-1"]][class_2_chicken_list[["0-1"]],5]
Table_S7[names(class_numbers_list[["0-1"]][class_2_chicken_list[["0-1"]],6]),'Human_Spec_Num'] <- class_numbers_list[["0-1"]][class_2_chicken_list[["0-1"]],6]

Table_S7[names(class_numbers_list[["0-1"]][class_2_both_list[["0-1"]],7]),'Conserved_Num'] <- class_numbers_list[["0-1"]][class_2_both_list[["0-1"]],7]
Table_S7[names(class_numbers_list[["0-1"]][class_2_both_list[["0-1"]],8]),'Human_Spec_Num'] <- class_numbers_list[["0-1"]][class_2_both_list[["0-1"]],8]
Table_S7[names(class_numbers_list[["0-1"]][class_2_both_list[["0-1"]],9]),'chicken_Spec_Num'] <- class_numbers_list[["0-1"]][class_2_both_list[["0-1"]],9]

Table_S7[names(class_numbers_list[["0-1"]][class_3_list[["0-1"]],10]),'Human_Spec_Num'] <- class_numbers_list[["0-1"]][class_3_list[["0-1"]],10]
Table_S7[names(class_numbers_list[["0-1"]][class_3_list[["0-1"]],11]),'chicken_Spec_Num'] <- class_numbers_list[["0-1"]][class_3_list[["0-1"]],11]

Table_S7[names(class_numbers_list[["0-1"]][conserved_exp_list[["0-1"]],12]),'Conserved_Num'] <- class_numbers_list[["0-1"]][conserved_exp_list[["0-1"]],12]

rownames(Table_S7) <- rownames(htiss.anno.comp@assays$RNA@data)



#####

Table_S7 <- matrix(0,ngenes,5)
rownames(Table_S7) <- rownames(htiss.anno.comp@assays$RNA@data)
colnames(Table_S7) <- c('Evo_Type', 'Gene_Class', 'Conserved_Num', 'Human_Spec_Num', 'chicken_Spec_Num')

Table_S7[conserved_exp_list[["0-0"]],'Evo_Type'] <- 0
Table_S7[conserved_not_exp_list[["0-0"]],'Evo_Type'] <- NA
Table_S7[class_1_human_list[["0-0"]],'Evo_Type'] <- 1
Table_S7[class_1_chicken_list[["0-0"]],'Evo_Type'] <- 1
Table_S7[class_2_human_list[["0-0"]],'Evo_Type'] <- 2
Table_S7[class_2_chicken_list[["0-0"]],'Evo_Type'] <- 2
Table_S7[class_2_both_list[["0-0"]],'Evo_Type'] <- 2
Table_S7[class_3_list[["0-0"]],'Evo_Type'] <- 3




Table_S7[names(class_numbers_list[["0-0"]][class_1_human_list[["0-0"]],1]),'Human_Spec_Num'] <- class_numbers_list[["0-0"]][class_1_human_list[["0-0"]],1]
Table_S7[names(class_numbers_list[["0-0"]][class_1_chicken_list[["0-0"]],2]),'chicken_Spec_Num'] <- class_numbers_list[["0-0"]][class_1_chicken_list[["0-0"]],2]

Table_S7[names(class_numbers_list[["0-0"]][class_2_human_list[["0-0"]],3]),'Conserved_Num'] <- class_numbers_list[["0-0"]][class_2_human_list[["0-0"]],3]
Table_S7[names(class_numbers_list[["0-0"]][class_2_human_list[["0-0"]],4]),'Human_Spec_Num'] <- class_numbers_list[["0-0"]][class_2_human_list[["0-0"]],4]

Table_S7[names(class_numbers_list[["0-0"]][class_2_chicken_list[["0-0"]],5]),'Conserved_Num'] <- class_numbers_list[["0-0"]][class_2_chicken_list[["0-0"]],5]
Table_S7[names(class_numbers_list[["0-0"]][class_2_chicken_list[["0-0"]],6]),'Human_Spec_Num'] <- class_numbers_list[["0-0"]][class_2_chicken_list[["0-0"]],6]

Table_S7[names(class_numbers_list[["0-0"]][class_2_both_list[["0-0"]],7]),'Conserved_Num'] <- class_numbers_list[["0-0"]][class_2_both_list[["0-0"]],7]
Table_S7[names(class_numbers_list[["0-0"]][class_2_both_list[["0-0"]],8]),'Human_Spec_Num'] <- class_numbers_list[["0-0"]][class_2_both_list[["0-0"]],8]
Table_S7[names(class_numbers_list[["0-0"]][class_2_both_list[["0-0"]],9]),'chicken_Spec_Num'] <- class_numbers_list[["0-0"]][class_2_both_list[["0-0"]],9]

Table_S7[names(class_numbers_list[["0-0"]][class_3_list[["0-0"]],10]),'Human_Spec_Num'] <- class_numbers_list[["0-0"]][class_3_list[["0-0"]],10]
Table_S7[names(class_numbers_list[["0-0"]][class_3_list[["0-0"]],11]),'chicken_Spec_Num'] <- class_numbers_list[["0-0"]][class_3_list[["0-0"]],11]

Table_S7[names(class_numbers_list[["0-0"]][conserved_exp_list[["0-0"]],12]),'Conserved_Num'] <- class_numbers_list[["0-0"]][conserved_exp_list[["0-0"]],12]

rownames(Table_S7) <- rownames(htiss.anno.comp@assays$RNA@data)


levels(htiss.anno.comp)
celltypes<-levels(htiss.anno.comp)
ngenes <- length(rownames(htiss.anno.comp@assays$RNA@data))
cnum <- length(celltypes)

average_cells.comp <- matrix(0,ngenes,cnum)
pct_exp.comp <- matrix(0,ngenes,cnum)

colnames(average_cells.comp) <- celltypes
rownames(average_cells.comp) <- rownames(htiss.anno.comp@assays$RNA@data)

colnames(pct_exp.comp) <- celltypes
rownames(pct_exp.comp) <- rownames(htiss.anno.comp@assays$RNA@data)


for (i in 1:cnum) {
  cells <- WhichCells(htiss.anno.comp, ident = celltypes[i])
  average_cells.comp[,i] <- rowMeans(as.matrix(htiss.anno.comp@assays$RNA@data[,cells]))
  average_cells.comp[which(is.nan(average_cells.comp))] <- 0
  pct_exp.comp[,i] <- rowSums(htiss.anno.comp@assays$RNA@data[,cells] > 0) / length(cells)
}


is_on.comp <- list()


for (j in 0) {
  for (k in 1) {
    
    tmp <- matrix(0,ngenes,cnum)
    
    colnames(tmp) <- celltypes
    rownames(tmp) <- rownames(htiss.anno.comp@assays$RNA@data)
    
    for (i in 1:cnum) {
      
      median_cutoff <- median(average_cells.comp[,i]) + j * sd(average_cells.comp[,i])
      pct_cutoff <- median(pct_exp.comp[,i]) + k * sd(pct_exp.comp[,i])
      
      tmp[intersect(names(which(average_cells.comp[,i] > median_cutoff)), names(which(pct_exp.comp[,i] > pct_cutoff))),i] <- 1
      
    }
    
    colnames(tmp) <- celltypes
    rownames(tmp) <- rownames(htiss.anno.comp@assays$RNA@data)
    
    is_on.comp[[paste0(j,'-',k)]] <- tmp
    
  }
}



chicken.celltypes<-levels(mtiss.anno.comp)
ngenes <- length(rownames(mtiss.anno.comp@assays$RNA@data))
chicken.cnum<-length(chicken.celltypes)





chicken.average_cells.comp <- matrix(0,ngenes,chicken.cnum)
chicken.pct_exp.comp <- matrix(0,ngenes,chicken.cnum)

colnames(chicken.average_cells.comp) <- chicken.celltypes
rownames(chicken.average_cells.comp) <- rownames(mtiss.anno.comp@assays$RNA@data)

colnames(chicken.pct_exp.comp) <- chicken.celltypes
rownames(chicken.pct_exp.comp) <- rownames(mtiss.anno.comp@assays$RNA@data)

for (i in 1:chicken.cnum) {
  cells <- WhichCells(mtiss.anno.comp, ident = chicken.celltypes[i])
  chicken.average_cells.comp[,i] <- rowMeans(as.matrix(mtiss.anno.comp@assays$RNA@data[,cells]))
  chicken.average_cells.comp[which(is.nan(chicken.average_cells.comp))] <- 0
  chicken.pct_exp.comp[,i] <- rowSums(mtiss.anno.comp@assays$RNA@data[,cells] > 0) / length(cells)
  
}


##########################
chicken.is_on.comp <- list()

for (j in 0) {
  for (k in 1) {
    
    tmp <- matrix(0,ngenes,chicken.cnum)
    
    colnames(tmp) <- chicken.celltypes
    rownames(tmp) <- rownames(mtiss.anno.comp@assays$RNA@data)
    
    for (i in 1:chicken.cnum) {
      
      median_cutoff <- median(chicken.average_cells.comp[,i]) + j * sd(chicken.average_cells.comp[,i])
      pct_cutoff <- median(chicken.pct_exp.comp[,i]) + k * sd(chicken.pct_exp.comp[,i])
      
      tmp[intersect(names(which(chicken.average_cells.comp[,i] > median_cutoff)), names(which(chicken.pct_exp.comp[,i] > pct_cutoff))),i] <- 1
      
    }
    
    colnames(tmp) <- chicken.celltypes
    rownames(tmp) <- rownames(mtiss.anno.comp@assays$RNA@data)
    
    chicken.is_on.comp[[paste0(j,'-',k)]] <- tmp
    
  }
}



