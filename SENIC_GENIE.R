install.package("BiocManager")
BiocManager::install(c("AUCell", "RcisTarget","GENIE3","zoo", "mixtools", 
                       "rbokeh","DT", "NMF", "pheatmap", "R2HTML", "Rtsne","doMC", "doRNG","scRNAseq"))
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC")
BiocManager::install("Seurat")
#check
library(SCENIC)
packageVersion("SCENIC")


library(Seurat)
library(SeuratData)

pbmc3k = readRDS("./pbmc3k.test.seurat.Rds")

write.csv(t(as.matrix(pbmc3k@assays$RNA@counts)),file = "for.scenic.data.csv")

conda activate pyscenic 
cat >change.py
#change.py 

import os,sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("for.scenic.data.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sample.loom",x.X.transpose(),row_attrs,col_attrs);


cat >scenic.bash



#2.1 grn
pyscenic grn \
--num_workers 10 \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
$tfs #转录因子文件，1839个基因的名字列表

#2.2 cistarget
pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 20  \
--mask_dropouts

#2.3 AUCell
pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 10




loom <- open_loom('out_SCENIC.loom') 
regulonAucThresholds <- get_regulon_thresholds(loom)
TF_list <- as.character(regulonAucThresholds)
TF_list <- gsub("\\(\\+\\)", "", TF_list)


exprMatr <- read.csv("scenic.data.csv")


rownames(exprMatr) <- exprMatr$X
exprMatr <- exprMatr[,2:ncol(exprMatr)]
exprMatr <- t(exprMatr)
exprMatr <- as.matrix(exprMatr)

rownames(exprMatr) <- gsub("\\.", "-", rownames(exprMatr))

regulators <- TF_list

weightMat <- GENIE3(exprMatr, regulators=regulators,nCores=10, nTrees=50)
linkList <- getLinkList(weightMat)

save.image(file="B.RData")
