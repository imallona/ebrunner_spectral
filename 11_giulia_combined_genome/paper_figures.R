rm(list=ls())
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(scDblFinder)
library(BiocManager)
library(tidyr)
library(gridExtra)
library(stringr)
library(reshape2)
library(Matrix)
library(readr)
library(rtracklayer)

opts_chunk$set(fig.width = 5,
               fig.height = 5,
               cache = TRUE,
               include = TRUE,
               dev = "png",
               cache.lazy = FALSE,
               warning = TRUE,
               message = TRUE)

options(bitmapType='cairo')

ID <- 'descriptive analysis RoCK and ROI'
BASE <- file.path('/home', 'gmoro','RStudio','paper') # how to I simlink file????????????????
WD <- file.path(BASE, ID)
DATA_WTA <- file.path(BASE, 'RoCK_ROIWTA/RoCK_ROIWTA/Solo.out/Gene/filtered/') # is it ok to used filtered??????????????
DATA_alien <- file.path(BASE, 'RoCK_ROIfeaturecounts/')
DATA_TSO <- file.path(BASE,'RoCK_ROITSO/RoCK_ROITSO/Solo.out/Gene/filtered/')
DOWNSAMPLE <- FALSE

## setwd(WD)

NTHREADS <- 10

plan("multisession", workers = NTHREADS)
options(future.globals.maxSize= 30e9)

ac <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}



#WTA ANALYSIS 



#WTA generating starsolo count matrix without alien chromosomes (those are in starsolo)




# loading WTA count table 

counts_RR_WTA <- readMM(file.path(DATA_WTA,'matrix.mtx'))
genes_RR_WTA <- read.table(file.path(DATA_WTA,'features.tsv'), header = FALSE)
gene_ids_RR_WTA <- genes_RR_WTA$V2
cell_ids_RR_WTA <- read.table(file.path(DATA_WTA,'barcodes.tsv'), header = FALSE)$X1

rownames(counts_RR_WTA) <- gene_ids_RR_WTA
colnames(counts_RR_WTA) <- cell_ids_RR_WTA

# removing alien chromosomes from WTA --> should we do this??????????????????????????

rownames(RR_WTA)[1:111]
RR_WTA<-RR_WTA[112:length(rownames(RR_WTA)),]
head(rownames(RR_WTA))





#Adding information from featurecounts table to alien chromosomes 

RR_table<-read.table(file.path(DATA_alien,'alien_only'),header=TRUE)

RR_table$Geneid[7] <- "capture_sequence_double_egfp_1"

rownames(RR_table) <- RR_table$Geneid
RR_table <- RR_table[1:20, -c(1:6)] 

RR_tso_table <- RR_table[,grep('tso', colnames(RR_table))]
RR_wta_table <- RR_table[,grep('wta', colnames(RR_table))]

# adapting colnames 

colnames(RR_wta_table)

b <- data.frame(barcodes=character()) # to optimize????????????????
for (i in 1:length(colnames(RR_wta_table))){
  new<-paste(strsplit(colnames(RR_wta_table),'_')[[i]][5],'_',strsplit(colnames(RR_wta_table),'_')[[i]][6],'_',strsplit(colnames(RR_wta_table),'_')[[i]][7],sep="")
  print(new)
  b[i,1]<-new
}

colnames(RR_wta_table) <- b$barcodes

saveRDS(RR_wta_table,file=file.path(BASE,'RR_WTA_table_correct.RData'))
RR_wta_table<-readRDS(file = file.path(BASE,'RR_WTA_table_correct.RData'))

# merging the two datasets --> sort and then merge --> how to do this better???????

RR_wta_table<-data.matrix(RR_wta_table[,sort(colnames(RR_wta_table))])
RR_WTA<-RR_WTA[,sort(colnames(RR_WTA))]

RR_WTA_combined <- rbind(RR_wta_table,RR_WTA)

saveRDS(RR_WTA_combined,file=file.path(BASE,'RR_WTA_combined.RData'))
RR_WTA_combined  <-readRDS(file = file.path(BASE,'RR_WTA_combined.RData'))

# changing the gene names 

