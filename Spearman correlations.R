library(dplyr)
library(Seurat)
library(patchwork)
library(future)
library(corrplot)
plan("multicore", workers = 10)
options(future.globals.maxSize= 62914560000)
##read one2one gene list
one2one <- read.csv("chicken_human_one2one.csv",header = T)


chi <- readRDS("chicken.rds")
hum <- readRDS("human.rds")

hum.raw <- GetAssayData(hum, slot="counts")
chi.raw <- GetAssayData(chi, slot="counts")
chi_inter <- chi.raw[which (row.names(chi.raw) %in% one2one$Chicken),]
inter_location <- match (row.names(chi_inter), one2one$Chicken)
row.names(chi_inter)  <- one2one$Human[inter_location]

hum_inter <- hum.raw[which (row.names(hum.raw) %in% one2one$Human),]

dup_gene <- row.names(chi_inter)[which(duplicated(row.names(chi_inter)))]

if(length(dup_gene)>0){
  chi_inter <- chi_inter[-which(row.names(chi_inter) %in% dup_gene),]
}
dup_gene2 <- row.names(hum_inter)[which(duplicated(row.names(hum_inter)))]
if(length(dup_gene2)>0){
  hum_inter <- hum_inter[-which(row.names(hum_inter) %in% dup_gene2)]
}



chi_inter <- na.omit(chi_inter)
hum_inter <- na.omit(hum_inter)


chicken <- CreateSeuratObject(chi_inter,meta.data = chi@meta.data)
human <- CreateSeuratObject(hum_inter,meta.data = hum@meta.data)


chicken <- NormalizeData(chicken, normalization.method = "LogNormalize", scale.factor = 10000)
chicken <- FindVariableFeatures(chicken, selection.method = "vst", nfeatures = 10000)
chicken <- ScaleData(chicken)

human <- NormalizeData(human, normalization.method = "LogNormalize", scale.factor = 10000)
human <- FindVariableFeatures(human, selection.method = "vst", nfeatures = 10000)
human <- ScaleData(human)


chickenmatrix <- AverageExpression(chicken,group.by = "cell_type", assays = "RNA")
chickenmatrix <- as.data.frame(chickenmatrix)

humanmatrix <- AverageExpression(human,group.by = "cell_type", assays = "RNA")
humanmatrix <- as.data.frame(humanmatrix)


chiDE=names(tail(sort(apply(chickenmatrix,1,sd)),2500))  
humDE=names(tail(sort(apply(humanmatrix,1,sd)),2500))  



chi_SD <- chickenmatrix[which(row.names(chickenmatrix) %in% chiDE),] 
hum_SD <- humanmatrix[which(row.names(humanmatrix) %in% humDE),]

nDEchi=length(row.names(chi_SD))
nDEhum=length(row.names(hum_SD))

common_genes <- intersect(chiDE,humDE)
chi_SD <- chi_SD[rownames(chi_SD) %in% common_genes,]   
hum_SD <- hum_SD[rownames(hum_SD) %in% common_genes,]


avg=rowMeans(chi_SD) 
chi_SD=sweep(chi_SD,1,avg,"/")   
avg=rowMeans(hum_SD)
hum_SD=sweep(hum_SD,1,avg,"/")
rm(avg)


colnames(chi_SD) <- gsub("RNA.","Chicken_",colnames(chi_SD))
colnames(hum_SD) <- gsub("RNA.","Human_",colnames(hum_SD)) 
geTable=merge(hum_SD,chi_SD,by='row.names', all=F) 
rownames(geTable) =geTable$Row.names  
geTable=geTable[,2:ncol(geTable)]
#
Corr.Coeff.Table=cor(geTable,method='spearman')  


#####################Visualization
shuffled.cor.list=list()
pb<-txtProgressBar(1, 100, style=3)

for(i in 1:500){
  shuffled=apply(geTable[,1:ncol(hum_SD)],1,sample)
  shuffled2=apply(geTable[,(ncol(hum_SD)+1):ncol(geTable)],1,sample)
  shuffled=cbind(t(shuffled),t(shuffled2))
  shuffled.cor=cor(shuffled,method='spearman')
  shuffled.cor.list[[i]] =shuffled.cor
  rm(list=c('shuffled','shuffled2','shuffled.cor'))
  if((i%%35) ==0){
    setTxtProgressBar(pb, (i*100)/1000)}}

######################################

p.value.table=matrix(ncol=ncol(geTable), nrow=ncol(geTable))
rownames(p.value.table) =colnames(geTable)
colnames(p.value.table) =colnames(geTable)

shuffled.mean.table=matrix(ncol=ncol(geTable), nrow=ncol(geTable))
rownames(shuffled.mean.table) =colnames(geTable)
colnames(shuffled.mean.table) =colnames(geTable)

a=combn(1:ncol(geTable),2)
for(i in 1:ncol(a)){
  cor.scores=sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
  shuffled.mean.table[a[1,i],a[2,i]] =mean(cor.scores)
  shuffled.mean.table[a[2,i],a[1,i]] =mean(cor.scores)
  cor.scores=sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
  p.value=mean(abs(cor.scores)>=abs(Corr.Coeff.Table[a[1,i],a[2,i]]))
  p.value.table[a[1,i],a[2,i]] =p.value
  p.value.table[a[2,i],a[1,i]] =p.value
  rm(list=c('cor.scores','p.value'))
  setTxtProgressBar(pb, (i*100)/ncol(a))}
#
neg.log10.p=-log10(p.value.table)
#
a = combn(1:ncol(geTable),2)
marker.overlap.list=list()
for(i in 1:ncol(a)){
  datasubset=cbind(geTable[,a[1,i]],geTable[,a[2,i]])
  markers=rownames(geTable[datasubset[,1]>1.5&datasubset[,2]>1.5,])
  marker.overlap.list[[i]] =markers
  names(marker.overlap.list)[i] =paste(colnames(geTable)[a[1,i]], colnames(geTable)[a[2,i]],sep='_')
  rm(list=c('datasubset','markers'))
}

ExpressionTableSpecies1 <- humanmatrix
DEgenesSpecies1 <- humDE
ExpressionTableSpecies2 <- chickenmatrix
DEgenesSpecies2 <- chiDE
Species1 <- "human"
Species2 <- "chicken"
filename <- "Corrplotend"
Height <- 5
Width <- 7

list.to.return=list(Corr.Coeff.Table,shuffled.mean.table,p.value.table,neg.log10.p,common_genes,rownames(geTable),
                    length(common_genes),length(rownames(geTable)),nDEchi,nDEhum,geTable)


names(list.to.return) =c('corr.coeff','shuffled_correlation_score_means','p.value',
                         'negative_log10_p.value','DEgenes_intersect','DEgenes_in_analysis',
                         'nDEgenes_intersect','nDEgenes_in_analysis','nDEgenes_chi','nDEgenes_hum','scaled_table')

comp.intersect <- list.to.return
rm(list.to.return)
#1
comp_table.intersect <- t(comp.intersect[[1]]
                          [1:ncol(ExpressionTableSpecies1),(ncol(ExpressionTableSpecies1)+1):nrow(comp.intersect[[1]])])
p_table.intersect <- t(comp.intersect[[3]][1:ncol(ExpressionTableSpecies1),(ncol(ExpressionTableSpecies1)+1):nrow(comp.intersect[[1]])])
p_table.intersect <- (-p_table.intersect)


pdf(paste0(filename, "_human_chi_1151_nocycling_1.0.pdf"), height=11, width = 13)
col1 <- colorRampPalette(c("darkblue", "white","darkred"))
par(mfrow=c(2,2), oma=c(0,0,2,0)+0.5)
#
library(corrplot)
corrplot(comp_table.intersect, order="original",tl.pos="lt", 
         method="color", tl.col="black",cl.lim=c(min(comp_table.intersect),max(comp_table.intersect)), 
         is.corr=F,tl.cex=0.7, sig.level=(-0.05),insig="pch", 
         pch=19, pch.cex=0.8,pch.col="black",p.mat=p_table.intersect, 
         col=col1(200), main= paste("spearsman",",",comp.intersect[[8]], "genes", sep=" "),
         cex.main=0.8, cl.align.text="l",mar = c(0, 0, 2, 0))  

#
mtext(paste(length(DEgenesSpecies1), Species1, "genes", "and", length(DEgenesSpecies2),  Species2,"genes,","all genes", sep=" "),outer=T)



dev.off()







