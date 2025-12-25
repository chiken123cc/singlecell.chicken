
###############################################metaneighbor
salloc -N 1 -n 64 -p xahcnormal -J d181
salloc -N 1 -n 64 -p xahcnormal -J d182
salloc -N 1 -n 64 -p xahcnormal -J d183
conda activate metaneighbor

library(Seurat)
pbmc<-readRDS('integrated_umap_6_3000_anno_final_4.0_umap.rds')


chi= pbmc[,pbmc@meta.data$species %in% c("chicken")]
hum= pbmc[,pbmc@meta.data$species %in% c("human")]
tur= pbmc[,pbmc@meta.data$species %in% c("turtle")]

dim(chi)
dim(hum)
dim(tur)

DefaultAssay(chi) <- "RNA"
DefaultAssay(hum) <- "RNA"
DefaultAssay(tur) <- "RNA"


Idents(chi)<-chi$general_celltype_1
Idents(hum)<-hum$general_celltype_1
Idents(tur)<-tur$general_celltype_1

table(hum$general_celltype_1)

chi<- subset(chi,downsample = 2000)
hum<- subset(hum,downsample = 2000)
tur<- subset(tur,downsample = 2000)

library(magrittr)
library(dplyr)

rm(pbmc)
rm(hum)
rm(tur)
gc()

chi<- PseudoCell(chi, "RNA","data","general_celltype_1",20)
saveRDS(chi,"chi.rds")

hum<- PseudoCell(hum, "RNA","data","general_celltype_1",20)
saveRDS(hum,"hum.rds")

tur<- PseudoCell(tur, "RNA","data","general_celltype_1",20)
saveRDS(tur,"tur.rds")

c <- chi@meta.data
c <- data.frame(Sample_ID = unlist(paste0("chi" , rownames(c))),
                Study_ID = rep("chi_tr",times = nrow(c)),
                Celltype = unlist(paste0("chi_" , c$orig.ident)))

h <- hum@meta.data
h <- data.frame(Sample_ID = unlist(paste0("hum" , rownames(h))),
                Study_ID = rep("hum_tr",times = nrow(h)),
                Celltype = unlist(paste0("hum_" , h$orig.ident)))

t <- tur@meta.data
t <- data.frame(Sample_ID = unlist(paste0("tur" , rownames(t))),
                Study_ID = rep("tur_tr",times = nrow(t)),
                Celltype = unlist(paste0("tur_" , t$orig.ident)))
P<-rbind(c,h,t)
colnames(P)<-c("Sample_ID","Study_ID","Celltype")


chid<-as.data.frame(chi@assays[["RNA"]]@data)
humd<-as.data.frame(hum@assays[["RNA"]]@data)
turd<-as.data.frame(tur@assays[["RNA"]]@data)

chidt<-chid
cn<- colnames(chidt)
colnames(chidt) <- paste0("chi",cn)

humdt<-humd
hn<- colnames(humdt)
colnames(humdt) <- paste0("hum",hn)

turdt<-turd
tn<- colnames(turdt)
colnames(turdt) <- paste0("tur",tn)


df1 <- merge(chidt, humdt, by = "row.names", all = FALSE)

rownames(df1)<-df1$Row.names
df1 <- df1[, -1]

df1 <- merge(df1, turdt, by = "row.names", all = FALSE)

rownames(df1)<-df1$Row.names
df1 <- df1[, -1]

####


data<-df1[,as.character(P$Sample_ID)]

celltypes <-unique(as.character(P$Celltype))

var.genes=get_variable_genes(data,P)

length(var.genes)

# run_MetaNeighbor_US
celltype.NV=run_MetaNeighbor_US(var.genes,data,celltypes,P)
write.table(celltype.NV,file="celltype.NV.out",sep="\t",quote=F)


rowcol = as.data.frame(rownames(celltype.NV))
species = c()
for (i in rowcol$`rownames(celltype.NV)`){
  species = append(species, c(strsplit(i, "_",  fixed = TRUE))[[1]][1])
}
rowcol$species = species
ann_colors=list(species=c(chi="#90719f",hum="#f0ac47",tur="#779043"))

rowcol1<-rowcol
row.names(rowcol1)<-rowcol1$`rownames(celltype.NV)`
rowcol1[, 1] <- NULL
library(pheatmap)
pheatmap(celltype.NV)
pheatmap(celltype.NV,annotation_col=rowcol1,annotation_row=rowcol1,annotation_colors=ann_colors)

###################

GatherData<- function(object,
                      assay,
                      slot_use,
                      ...) {
  
  assay <- assay %||% "RNA"
  slot_use <- slot_use %||% "data"
  obj_data <- GetAssayData(
    object = object,
    assay = assay,
    slot = slot_use
  ) %>%
    as.matrix()
  return(obj_data)
}



PseudoCell  <- function(object,
                        assay_use = NULL,
                        slot_use = NULL,
                        cluster_use =NULL,
                        pseudocell.size  =NULL){
  message("tips: 
  Cluster_use : one col in metadata
  pseudocell.size : how many cell will be pseudo")
  
  Inter<- GatherData(object = object,
                     assay = assay_use,
                     slot_use = slot_use) 
  Inter[Inter<0]=0
  idd<-object@meta.data
  Inter.id<-cbind(rownames(idd),as.vector(idd[,cluster_use]))
  
  rownames(Inter.id)<-rownames(idd)
  colnames(Inter.id)<-c("CellID","Celltype")
  
  Inter.id<-as.data.frame(Inter.id)
  Inter1<-Inter[,Inter.id$CellID]
  Inter<-as.matrix(Inter1)
  pseudocell.size = pseudocell.size ## 10 test
  new_ids_list = list()
  Inter.id$Celltype <- as.factor(Inter.id$Celltype)
  for (i in 1:length(levels(Inter.id$Celltype))) {
    cluster_id = levels(Inter.id$Celltype)[i]
    cluster_cells <- rownames(Inter.id[Inter.id$Celltype == cluster_id,])
    cluster_size <- length(cluster_cells)       
    pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
    pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
    names(pseudo_ids) <- sample(cluster_cells)  
    new_ids_list[[i]] <- pseudo_ids     
  }
  
  new_ids <- unlist(new_ids_list)
  new_ids <- as.data.frame(new_ids)
  new_ids_length <- table(new_ids)
  
  new_colnames <- rownames(new_ids)  ###add
  all.data<-Inter[,as.character(new_colnames)] ###add
  all.data <- t(all.data)###add
  
  new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                      list(name=new_ids[,1]),FUN=mean)
  rownames(new.data)<-new.data$name
  new.data<-new.data[,-1]
  
  new_ids_length<-as.matrix(new_ids_length)##
  short<-which(new_ids_length< pseudocell.size -1 )##
  new_good_ids<-as.matrix(new_ids_length[-short,])##
  result<-t(new.data)[,rownames(new_good_ids)]
  rownames(result)<-rownames(Inter)
  
  newobject <- CreateSeuratObject(result)
  newobject@misc$idtrans <- new_ids
  newobject@commands$PseudoCell <- LogSeuratCommand(newobject, return.command = TRUE)
  return(newobject)
}



run_MetaNeighbor_US<-function(vargenes, data, celltypes, pheno){
  
  cell.labels=matrix(0,ncol=length(celltypes),nrow=dim(pheno)[1])
  rownames(cell.labels)=colnames(data)
  colnames(cell.labels)=celltypes
  for(i in 1:length(celltypes)){
    type=celltypes[i]
    m<-match(pheno$Celltype,type)
    cell.labels[!is.na(m),i]=1
  }
  
  m<-match(rownames(data),vargenes)
  cor.dat=cor(data[!is.na(m),],method="s")
  rank.dat=cor.dat*0
  rank.dat[]=rank(cor.dat,ties.method="average",na.last = "keep")
  rank.dat[is.na(rank.dat)]=0
  rank.dat=rank.dat/max(rank.dat)
  sumin    =  (rank.dat) %*% cell.labels
  sumall   = matrix(apply(rank.dat,2,sum), ncol = dim(sumin)[2], nrow=dim(sumin)[1])
  predicts = sumin/sumall
  
  cell.NV=matrix(0,ncol=length(celltypes),nrow=length(celltypes))
  colnames(cell.NV)=colnames(cell.labels)
  rownames(cell.NV)=colnames(cell.labels)
  
  for(i in 1:dim(cell.labels)[2]){
    predicts.temp=predicts
    m<-match(pheno$Celltype,colnames(cell.labels)[i])
    study=unique(pheno[!is.na(m),"Study_ID"])
    m<-match(pheno$Study_ID,study)
    pheno2=pheno[!is.na(m),]
    predicts.temp=predicts.temp[!is.na(m),]
    predicts.temp=apply(abs(predicts.temp), 2, rank,na.last="keep",ties.method="average")
    filter=matrix(0,ncol=length(celltypes),nrow=dim(pheno2)[1])
    m<-match(pheno2$Celltype,colnames(cell.labels)[i])
    filter[!is.na(m),1:length(celltypes)]=1
    negatives = which(filter == 0, arr.ind=T)
    positives = which(filter == 1, arr.ind=T)
    predicts.temp[negatives] <- 0
    np = colSums(filter,na.rm=T)
    nn = apply(filter,2,function(x) sum(x==0,na.rm=T))
    p =  apply(predicts.temp,2,sum,na.rm=T)
    cell.NV[i,]= (p/np - (np+1)/2)/nn
  }
  
  cell.NV=(cell.NV+t(cell.NV))/2
  return(cell.NV)
  
}

get_variable_genes<-function(data, pheno) {
  var.genes1=vector("list")
  experiment=unique(pheno$Study_ID)
  j=1
  for(exp in experiment){
    dat.sub=data[,pheno$Study_ID==exp]
    genes.list=vector("list")
    med.dat=apply(dat.sub,1,median)
    var.dat=apply(dat.sub,1,var)
    quant.med=unique(quantile(med.dat,prob=seq(0,1,length=11),type=5))
    genes.list=vector("list",length=length(quant.med))
    for(i in 1:length(quant.med)){
      if(i==1){
        filt1=med.dat<=quant.med[i]
        var.temp=var.dat[filt1]
        quant.var=quantile(var.temp,na.rm=T)
        filt2=var.temp>quant.var[4]###### total is 4;TF is3
        genes.list[[i]]=names(var.temp)[filt2]
      }
      else {
        filt1=med.dat<=quant.med[i]&med.dat>quant.med[i-1]
        var.temp=var.dat[filt1]
        quant.var=quantile(var.temp,na.rm=T)
        filt2=var.temp>quant.var[4]######
        genes.list[[i]]=names(var.temp)[filt2]
      }
    }
    temp=length(genes.list)
    var.genes1[[j]]=unlist(genes.list[1:temp-1])
    j=j+1
  }
  var.genes=Reduce(intersect, var.genes1)
  return(var.genes)
}


get_top_hits <- function(cell.NV, pheno, threshold=0.95, filename) {
  
  type_by_study=table(pheno[,c("Celltype","Study_ID")])
  m<-match(rownames(cell.NV),rownames(type_by_study))
  f.a=!is.na(m)
  f.b=m[f.a]
  cell.NV=cell.NV[f.a,f.a]
  type_by_study=type_by_study[f.b,]
  
  for(i in 1:dim(type_by_study)[2]){
    filt=type_by_study[,i]!=0
    cell.NV[filt,filt]=0
  }
  
  diag(cell.NV)=0
  temp=vector()
  for(i in 1:dim(cell.NV)[1]){
    temp=c(temp,which.max(cell.NV[i,]))
  }
  temp=cbind(rownames(cell.NV),temp)
  for(i in 1:dim(cell.NV)[1]){
    temp[i,2]=cell.NV[i,as.numeric(temp[i,2])]
  }
  
  recip=temp[duplicated(temp[,2]),]
  filt=as.numeric(temp[,2])>=threshold
  recip=rbind(recip,temp[filt,])
  recip=cbind(recip,c(rep("Reciprocal_top_hit",each=dim(recip)[1]-sum(filt)),rep(paste("Above",threshold,sep="_"),each=sum(filt))))
  recip=recip[!duplicated(recip[,2]),]
  
  recip2=cbind(rownames(recip),recip[,1:3])
  colnames(recip2)=c("Celltype_1","Celltype_2","Mean_AUROC","Match_type")
  rownames(recip2)=NULL
  recip=recip2[order(recip2[,3],decreasing=T),]
  recip2=as.data.frame(recip)
  recip2[,3]=round(as.numeric(as.character(recip2[,3])),2)
  write.table(recip,file=filename,sep="\t",quote=F)
  return(recip2)
}


get_meta_matrix <- function(rds, g2mg, species) {
  
  ## 
  counts <- GetAssayData(rds, slot="counts", assay="RNA")
  origin_genes = rownames(rds@assays$RNA)
  orth_genes = intersect(g2mg[, 1], origin_genes)
  counts.sub <- counts[orth_genes, ]
  
  ## 
  g2mg_hash = hash(keys = c(g2mg[, 1]), values = c(g2mg[, 2]))
  mega_genes = c()
  for (c in orth_genes){
    mega_genes = append(mega_genes, g2mg_hash[[c]])
  }
  row.names(counts.sub) = mega_genes
  # write.csv(as.matrix(counts.sub), file=paste(species, "_meta_gene.csv", sep=""), , quote=F)
  ##
  counts.sub = aggregate(counts.sub, by=list(mega_genes), FUN=sum)
  rownames(counts.sub) <- counts.sub$`Group.1`
  counts.sub = subset(counts.sub, select = -c(`Group.1`))
  new_seurat_object <- CreateSeuratObject(counts=counts.sub, meta.data = rds@meta.data)
  saveRDS(new_seurat_object, file=paste(species, "_meta_gene.rds", sep=""))
}

