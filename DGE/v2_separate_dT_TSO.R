### clearing environment 

rm(list=ls())

### loading packages 

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(scDblFinder)
library(BiocManager)
library(tidyr)
library(gridExtra)
library(stringr)
library(gplots)
library(RColorBrewer)
library(knitr)
library(scater)
library(data.table)
library(SingleCellExperiment)
library(scran)
library(dplyr)
library(DT)
library(pheatmap)
library(Cairo)
library(RColorBrewer)
library(future)
library(BiocParallel)
library(scDblFinder)
library(scSorter)
library(intrinsicDimension)
library(UpSetR)

# library(biomaRt)

### loading SB data

dt_100_100 <-read.csv("/home/gmoro/RStudio/SPECTRAL/sixth_scRNAseq_separate/dT_100_100_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")

colnames(dt_100_100)

which(colnames(dt_100_100)=="egfp") #18670

t_dt_100_100 <-t(dt_100_100)

length(colnames(t_dt_100_100)) # number of cells: 19595



dt_unmod <-read.csv("/home/gmoro/RStudio/SPECTRAL/sixth_scRNAseq_separate/dT_unmod_no_cap_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")

colnames(dt_unmod)

which(colnames(dt_unmod)=="egfp") #19072

t_dt_unmod <-t(dt_unmod)

length(colnames(t_dt_unmod)) # number of cells: 20268


### Venn diagram genes

genes_dt_100_100 <- rownames(t_dt_100_100)
genes_unmod <- rownames(t_dt_unmod)

length(genes_dt_100_100) # 23163
length(genes_unmod) # 23190

length(intersect(genes_dt_100_100,genes_unmod)) # 21469

human_unmod_comb<-list(genes_dt_100_100=genes_dt_100_100,genes_unmod=genes_unmod)

venn(human_unmod_comb)

### combining Seurat objects

S_dt_100_100 <- CreateSeuratObject(counts = t_dt_100_100, min.cells = 3, min.features = 200, project="mod")
S_dt_unmod <- CreateSeuratObject(counts = t_dt_unmod, min.cells = 3, min.features = 200, project="unmod")

combined <- merge(S_dt_100_100, y = S_dt_unmod, add.cell.ids = c("100_100", "unmod"), project = "unmod_vs_mod")
combined


combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt")

grep( "^mt", rownames(combined@assays$RNA), value = T)

VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,cols="violet") 

### UMAP

combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(combined), 10)

all.genes <- rownames(combined)

combined <- ScaleData(combined, features = all.genes)

combined <- RunPCA(combined, features = VariableFeatures(object = combined))

combined <- FindNeighbors(combined, dims = 1:10)
combined <- FindClusters(combined, resolution = 0.5)

combined <- RunUMAP(combined, dims = 1:10)

DimPlot(combined, reduction = "umap",group.by = "orig.ident")

### differential expression analysis

Idents(object = combined) <- "orig.ident"
markers <- FindMarkers(object = combined, ident.1="unmod",ident.2="mod")

length(rownames(markers)) # 2304

nrow(combined) #19715 --> total number of genes

table(combined$orig.ident)















### same thing with unmod and all mods


### loading SB data

dt_10_100 <-read.csv("/home/gmoro/RStudio/SPECTRAL/sixth_scRNAseq_separate/dT_10_100_no_cap_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")

colnames(dt_10_100)

which(colnames(dt_10_100)=="egfp") #18670

t_dt_10_100 <-t(dt_10_100)

length(colnames(t_dt_10_100)) # number of cells: 19595



dt_100_10 <-read.csv("/home/gmoro/RStudio/SPECTRAL/sixth_scRNAseq_separate/dT_100_10_no_cap_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")

colnames(dt_100_10)

which(colnames(dt_100_10)=="egfp") #18670

t_dt_100_10 <-t(dt_100_10)

length(colnames(t_dt_100_10)) # number of cells: 19595


### Venn diagram genes

genes_dt_100_100 <- rownames(t_dt_100_100)
genes_unmod <- rownames(t_dt_unmod)
genes_dt_100_10 <- rownames(t_dt_100_10)
genes_dt_10_100 <- rownames(t_dt_10_100)

cells_dt_100_100 <- colnames(t_dt_100_100)
cells_unmod <- colnames(t_dt_unmod)
cells_dt_100_10 <- colnames(t_dt_100_10)
cells_dt_10_100 <- colnames(t_dt_10_100)


length(genes_dt_100_100) # 23163
length(genes_unmod) # 23190
length(genes_dt_100_10) # 22885
length(genes_dt_10_100) # 21475

length(cells_dt_100_100) # 19738
length(cells_unmod) # 19894
length(cells_dt_100_10) # 10730
length(cells_dt_10_100) # 9882 

#length(intersect(genes_dt_100_100,genes_unmod)) # 21469

human_unmod_comb<-list(genes_dt_100_100=genes_dt_100_100,genes_unmod=genes_unmod,genes_dt_100_10=genes_dt_100_10,genes_dt_10_100=genes_dt_10_100)

venn(human_unmod_comb)

### combining Seurat objects

S_dt_100_100 <- CreateSeuratObject(counts = t_dt_100_100, min.cells = 3, min.features = 200, project="100_100_mod")
S_dt_unmod <- CreateSeuratObject(counts = t_dt_unmod, min.cells = 3, min.features = 200, project="unmod")
S_dt_100_10 <- CreateSeuratObject(counts = t_dt_100_10, min.cells = 3, min.features = 200, project="100_10_mod")
S_dt_10_100 <- CreateSeuratObject(counts = t_dt_10_100, min.cells = 3, min.features = 200, project="10_100_mod")

combined <- merge(S_dt_100_100, y = c(S_dt_unmod,S_dt_100_10,S_dt_10_100), add.cell.ids = c("100_100", "unmod","100_10","10_100"), project = "unmod_vs_mod")
combined

table(combined$orig.ident)

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt")

grep( "^mt", rownames(combined@assays$RNA), value = T)

VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

### UMAP

combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(combined), 10)

all.genes <- rownames(combined)

combined <- ScaleData(combined, features = all.genes)

combined <- RunPCA(combined, features = VariableFeatures(object = combined))

combined <- FindNeighbors(combined, dims = 1:10)
combined <- FindClusters(combined, resolution = 0.5)

combined <- RunUMAP(combined, dims = 1:10)

DimPlot(combined, reduction = "umap",group.by = "orig.ident")

### differential expression analysis

Idents(object = combined) <- "orig.ident"
markers_unmod_100_100 <- FindMarkers(object = combined, ident.1="unmod",ident.2="100_100_mod")
markers_unmod_100_10 <- FindMarkers(object = combined, ident.1="unmod",ident.2="100_10_mod")
markers_unmod_10_100 <- FindMarkers(object = combined, ident.1="unmod",ident.2="10_100_mod")

markers_100_100_10_100 <- FindMarkers(object = combined, ident.1="100_100_mod",ident.2="10_100_mod")
markers_100_100_100_10 <- FindMarkers(object = combined, ident.1="100_100_mod",ident.2="100_10_mod")
markers_10_100_100_10 <- FindMarkers(object = combined, ident.1="10_100_mod",ident.2="100_10_mod")


length(rownames(markers_10_100_100_10)) # 2304

nrow(combined) #20170 --> total number of genes

table(combined$orig.ident)

m_unmod_100_10 <-rownames(markers_unmod_100_10)
m_unmod_100_100 <-rownames(markers_unmod_100_100)
m_unmod_10_100 <-rownames(markers_unmod_10_100)

m_unmod_vs_mod <- intersect(intersect(m_unmod_100_10,m_unmod_100_100),m_unmod_10_100) 

length(m_unmod_vs_mod)

sort(m_unmod_vs_mod)

### venn diagram

human_unmod_comb<-list(unmod_100_10=m_unmod_100_10, unmod_100_100=m_unmod_100_100, unmod_10_100=m_unmod_10_100)

venn(human_unmod_comb)

### checking which genes actually expressed in samples or which are genes which are only expressed in one sample but lowly expressed

length(m_unmod_100_100)

length(intersect(m_unmod_100_100,rownames(t_dt_100_100))) #all genes in mod sample 
length(intersect(m_unmod_100_100,rownames(t_dt_unmod))) # all genes in unmod sample

length(m_unmod_100_10)

length(intersect(m_unmod_100_10,rownames(t_dt_100_10))) #all genes in mod sample 
length(intersect(m_unmod_100_10,rownames(t_dt_unmod))) # all genes in unmod sample

length(m_unmod_10_100)

length(intersect(m_unmod_10_100,rownames(t_dt_10_100))) #all genes in mod sample 
length(intersect(m_unmod_10_100,rownames(t_dt_unmod))) # all genes in unmod sample


### finding which genes have expression level more than 5 in a given dataset and refining list 

# finding for genes number of cells with gene expression more than 3

### for t_dt_100_100 --> larger than 3 and more than 50% cells so 9869

exp_100_100<-c()

for (i in m_unmod_vs_mod){
  #print(i)
  if (length(which(t_dt_100_100[which(rownames(t_dt_100_100)==i),]>3))>9869){
    print(i)
    exp_100_100<-c(exp_100_100,i)
  }
}

length(exp_100_100)

### for t_dt_100_10 --> larger than 3 and more than 50% cells so 5365

exp_100_10<-c()

for (i in m_unmod_vs_mod){
  #print(i)
  if (length(which(t_dt_100_10[which(rownames(t_dt_100_10)==i),]>3))>5365){
    print(i)
    exp_100_10<-c(exp_100_10,i)
  }
}

length(exp_100_10)

### for t_dt_10_100 --> larger than 3 and more than 50% cells so 4941

exp_10_100<-c()

for (i in m_unmod_vs_mod){
  #print(i)
  if (length(which(t_dt_10_100[which(rownames(t_dt_10_100)==i),]>3))>4941){
    print(i)
    exp_10_100<-c(exp_10_100,i)
  }
}

length(exp_10_100)

### for unmod --> larger than 3 and more than 50% cells so 9947

exp_unmod <-c()

for (i in m_unmod_vs_mod){
  #print(i)
  if (length(which(t_dt_unmod[which(rownames(t_dt_unmod)==i),]>3))>9947){
    print(i)
    exp_unmod<-c(exp_unmod,i)
  }
}

length(exp_unmod)









### rerun FindMarkers with highly expressed markers 

### for unmod --> larger than 3 and more than 50% cells so 9947

high_unmod <- c()

for (i in rownames(t_dt_unmod)){
  #print(i)
  if (length(which(t_dt_unmod[which(rownames(t_dt_unmod)==i),]>3))>9947){
    #print(i)
    high_unmod<-c(high_unmod,i)
  }
}

length(high_unmod)

m_high_unmod <- t_dt_unmod[high_unmod,]

### for t_dt_100_100 --> larger than 3 and more than 50% cells so 9869

high_100_100 <- c()

for (i in rownames(t_dt_100_100)){
  #print(i)
  if (length(which(t_dt_100_100[which(rownames(t_dt_100_100)==i),]>3))>9869){
    #print(i)
    high_100_100<-c(high_100_100,i)
  }
}

#high_100_100
length(high_100_100)

m_high_100_100 <- t_dt_100_100[high_100_100,]

### for t_dt_100_10 --> larger than 3 and more than 50% cells so 5365


high_100_10 <- c()

for (i in rownames(t_dt_100_10)){
  print(i)
  if (length(which(t_dt_100_10[which(rownames(t_dt_100_10)==i),]>3))>5365){
    #print(i)
    high_100_10<-c(high_100_10,i)
  }
}

#high_100_10
length(high_100_10)

m_high_100_10 <- t_dt_100_10[high_100_10,]


### for t_dt_10_100 --> larger than 3 and more than 50% cells so 4941


high_10_100 <- c()

for (i in rownames(t_dt_10_100)){
  #print(i)
  if (length(which(t_dt_10_100[which(rownames(t_dt_10_100)==i),]>3))>4941){
    #print(i)
    high_10_100<-c(high_10_100,i)
  }
}

#high_10_100
length(high_10_100)

m_high_10_100 <- t_dt_10_100[high_10_100,]

### example for a few highly expressed genes










### QC with codes from Iza

# converting int class Single Cell Experiment

combined_sc <- as.SingleCellExperiment(combined)

combined_sc <- addPerCellQC(combined_sc, 
                            subsets=list(Mito=grep("mt-", rownames(combined_sc))))


# plots 

plotColData(combined_sc, x = "sum", y="detected", colour_by="orig.ident") 

plotHighestExprs(combined_sc, exprs_values = "counts")

combined_sc <- logNormCounts(combined_sc)

vars <- getVarianceExplained(combined_sc,
                             variables="orig.ident")

head(vars)

plotExplanatoryVariables(vars)

# additional QC plots

libsize.drop <- isOutlier(combined_sc$total, nmads = 2, type = "both", log = TRUE,
                          batch = combined_sc$orig.ident)

feature.drop <- isOutlier(combined_sc$detected, nmads = 2, type = "both", log = TRUE,
                          batch = combined_sc$orig.ident)

mito.drop <- isOutlier(combined_sc$subsets_Mito_percent, nmads = 2, type = "higher",
                       batch = combined_sc$orig.ident, log = FALSE)


libsize.drop.manual <- combined_sc$total < 5000

par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))

curr <- subset(combined_sc, ,orig.ident == "100_100_mod")

hist(curr$total/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")








